#include "em.h"
#include "classify.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

em_config_t em_config_default(void) {
    return (em_config_t){
        .max_iter = 200,
        .conv_threshold = 1e-6,
        .n_restarts = 3,
        .alpha = 0.5,          /* Dirichlet sparsity prior */
        .d_mu = 0.0,
        .d_sigma = 0.5,       /* LogNormal prior on DNA yield */
        .b_mu = 0.0,
        .b_sigma = 0.5,       /* LogNormal prior on PCR bias */
        .estimate_degradation = 0,
        .min_containment = 0.1,
        .seed = 42,
        .mito_copy_numbers = NULL,
        .prune_threshold = 0.0,
        .use_advanced_ci = 0,
        .use_brent_lambda = 0,
        .use_full_lrt = 0,
    };
}

/* --- Internal EM state --- */
typedef struct {
    double *w;     /* [S] weight fractions */
    double *d;     /* [S] DNA yield factors */
    double *b;     /* [S*M] PCR bias (row-major by species) */
    double lambda;
    int S, M;
    const int *amp_lens;         /* [S*M] amplicon lengths (borrowed, not owned) */
    int estimate_degradation;    /* flag: include exp(-lambda*L) term */
} em_params_t;

/* Forward declaration for Brent's method (used in m_step) */
static double brent_lambda_update(const em_params_t *p, const em_read_t *reads,
                                   int n_reads, double **gamma);

static em_params_t *params_alloc(int S, int M) {
    em_params_t *p = (em_params_t *)hs_calloc(1, sizeof(em_params_t));
    p->S = S; p->M = M;
    p->w = (double *)hs_calloc((size_t)S, sizeof(double));
    p->d = (double *)hs_calloc((size_t)S, sizeof(double));
    p->b = (double *)hs_calloc((size_t)(S * M), sizeof(double));
    return p;
}

static void params_free(em_params_t *p) {
    if (!p) return;
    free(p->w); free(p->d); free(p->b); free(p);
}

static void params_init_random(em_params_t *p, hs_rng_t *rng, const em_config_t *cfg) {
    int S = p->S, M = p->M;
    /* Random Dirichlet initialization for w */
    double *alpha = (double *)hs_malloc((size_t)S * sizeof(double));
    for (int s = 0; s < S; s++) alpha[s] = cfg->alpha + 0.5;
    hs_rng_dirichlet(rng, alpha, p->w, S);
    free(alpha);

    /* Initialize d near 1 */
    for (int s = 0; s < S; s++)
        p->d[s] = hs_rng_lognormal(rng, cfg->d_mu, 0.1);

    /* Initialize b near 1 */
    for (int s = 0; s < S; s++)
        for (int m = 0; m < M; m++)
            p->b[s * M + m] = hs_rng_lognormal(rng, cfg->b_mu, 0.1);

    p->lambda = 0.001;
}

static void params_init_uniform(em_params_t *p) {
    int S = p->S, M = p->M;
    for (int s = 0; s < S; s++) p->w[s] = 1.0 / S;
    for (int s = 0; s < S; s++) p->d[s] = 1.0;
    for (int s = 0; s < S; s++)
        for (int m = 0; m < M; m++) p->b[s * M + m] = 1.0;
    p->lambda = 0.001;
}

/* --- E-step: compute responsibilities --- */
static double e_step(const em_params_t *p, const em_read_t *reads, int n_reads,
                     double **gamma) {
    int M = p->M;
    double ll = 0.0;

    for (int r = 0; r < n_reads; r++) {
        int nc = reads[r].n_candidates;
        int m = reads[r].marker_idx;
        if (m < 0) m = 0; /* fallback */

        double max_log = -INFINITY;
        double *log_g = (double *)hs_malloc((size_t)nc * sizeof(double));

        for (int j = 0; j < nc; j++) {
            int s = reads[r].species_indices[j];
            double c = reads[r].containments[j];
            if (c < 1e-300) c = 1e-300;
            log_g[j] = log(p->w[s]) + log(p->d[s]) + log(p->b[s * M + m]) + log(c);
            if (p->estimate_degradation && p->amp_lens) {
                int L = p->amp_lens[s * M + m];
                if (L > 0) log_g[j] -= p->lambda * (double)L;
            }
            if (log_g[j] > max_log) max_log = log_g[j];
        }

        /* logsumexp normalization */
        double sum = 0.0;
        for (int j = 0; j < nc; j++)
            sum += exp(log_g[j] - max_log);
        double log_norm = max_log + log(sum);

        for (int j = 0; j < nc; j++)
            gamma[r][j] = exp(log_g[j] - log_norm);

        ll += log_norm;
        free(log_g);
    }

    return ll;
}

/* --- M-step: update parameters --- */
static void m_step(em_params_t *p, const em_read_t *reads, int n_reads,
                   double **gamma, const em_config_t *cfg, int single_marker) {
    int S = p->S, M = p->M;
    double *eff_counts = (double *)hs_calloc((size_t)S, sizeof(double));

    /* Accumulate effective counts for w */
    for (int r = 0; r < n_reads; r++) {
        for (int j = 0; j < reads[r].n_candidates; j++) {
            int s = reads[r].species_indices[j];
            eff_counts[s] += gamma[r][j];
        }
    }

    /* Update w with Dirichlet MAP.
     * In single-marker mode, d is fixed to mito CN priors, so we divide out
     * the DNA yield effect: w_s ∝ eff_counts_s / d_s. This corrects for
     * species that produce more reads due to higher mito copy number. */
    double w_sum = 0.0;
    for (int s = 0; s < S; s++) {
        double adj = single_marker && p->d[s] > 1e-10 ? p->d[s] : 1.0;
        p->w[s] = (eff_counts[s] / adj) + cfg->alpha - 1.0;
        if (p->w[s] < 1e-10) p->w[s] = 1e-10;
        w_sum += p->w[s];
    }
    for (int s = 0; s < S; s++) p->w[s] /= w_sum;

    /* In single-marker mode, d and b are fixed — skip their updates */
    if (single_marker) {
        free(eff_counts);
        return;
    }

    /* Update d: DNA yield factors
     * d_s is proportional to (effective reads for s) / (w_s * sum_m b_sm * reads_at_m)
     * With LogNormal MAP: d_s = exp( (N_s/sigma^2 + mu/sigma_prior^2) / (N_s/sigma^2 + 1/sigma_prior^2) )
     * Simplified: weighted ratio approach */
    double *d_num = (double *)hs_calloc((size_t)S, sizeof(double));
    double *d_den = (double *)hs_calloc((size_t)S, sizeof(double));
    for (int r = 0; r < n_reads; r++) {
        int m = reads[r].marker_idx;
        if (m < 0) m = 0;
        for (int j = 0; j < reads[r].n_candidates; j++) {
            int s = reads[r].species_indices[j];
            d_num[s] += gamma[r][j];
            d_den[s] += p->w[s] * p->b[s * M + m];
        }
    }
    for (int s = 0; s < S; s++) {
        if (d_den[s] > 1e-10) {
            double ratio = d_num[s] / d_den[s];
            /* LogNormal MAP regularization */
            double log_ratio = log(ratio > 1e-10 ? ratio : 1e-10);
            double prior_prec = 1.0 / (cfg->d_sigma * cfg->d_sigma);
            double data_prec = eff_counts[s] > 1.0 ? eff_counts[s] : 1.0;
            double post_mean = (data_prec * log_ratio + prior_prec * cfg->d_mu) /
                              (data_prec + prior_prec);
            p->d[s] = exp(post_mean);
        }
    }
    /* Normalize d: geometric mean = 1 (identifiability) */
    double log_d_sum = 0.0;
    int d_count = 0;
    for (int s = 0; s < S; s++) {
        if (p->d[s] > 1e-10) { log_d_sum += log(p->d[s]); d_count++; }
    }
    if (d_count > 0) {
        double geomean = exp(log_d_sum / d_count);
        for (int s = 0; s < S; s++) p->d[s] /= geomean;
    }
    free(d_num); free(d_den);

    /* Update b: PCR bias per species x marker
     * b_sm proportional to (reads assigned to s from marker m) / (d_s * w_s * total_at_m) */
    double *b_num = (double *)hs_calloc((size_t)(S * M), sizeof(double));
    double *marker_total = (double *)hs_calloc((size_t)M, sizeof(double));
    for (int r = 0; r < n_reads; r++) {
        int m = reads[r].marker_idx;
        if (m < 0) m = 0;
        for (int j = 0; j < reads[r].n_candidates; j++) {
            int s = reads[r].species_indices[j];
            b_num[s * M + m] += gamma[r][j];
        }
        marker_total[m] += 1.0;
    }
    for (int s = 0; s < S; s++) {
        for (int m = 0; m < M; m++) {
            double expected = p->w[s] * p->d[s] * marker_total[m];
            if (expected > 1e-10 && b_num[s * M + m] > 0) {
                double ratio = b_num[s * M + m] / expected;
                double log_ratio = log(ratio > 1e-10 ? ratio : 1e-10);
                double prior_prec = 1.0 / (cfg->b_sigma * cfg->b_sigma);
                double data_prec = b_num[s * M + m] > 1.0 ? b_num[s * M + m] : 1.0;
                double post_mean = (data_prec * log_ratio + prior_prec * cfg->b_mu) /
                                  (data_prec + prior_prec);
                p->b[s * M + m] = exp(post_mean);
            } else {
                /* Prior only */
                p->b[s * M + m] = exp(cfg->b_mu);
            }
        }
        /* Geometric mean normalization per species (identifiability) */
        double log_b_sum = 0.0;
        int b_count = 0;
        for (int m = 0; m < M; m++) {
            if (p->b[s * M + m] > 1e-10) {
                log_b_sum += log(p->b[s * M + m]);
                b_count++;
            }
        }
        if (b_count > 0) {
            double geomean = exp(log_b_sum / b_count);
            for (int m = 0; m < M; m++) p->b[s * M + m] /= geomean;
        }
    }
    free(b_num); free(marker_total);

    /* Update lambda: degradation rate */
    if (cfg->estimate_degradation && p->amp_lens) {
        if (cfg->use_brent_lambda) {
            /* Brent's method: maximize Q(lambda) exactly */
            p->lambda = brent_lambda_update(p, reads, n_reads, gamma);
        } else {
            /* Closed-form moment matching */
            double gamma_L_sum = 0.0;
            double gamma_sum = 0.0;
            for (int r = 0; r < n_reads; r++) {
                int m_idx = reads[r].marker_idx;
                if (m_idx < 0) m_idx = 0;
                for (int j = 0; j < reads[r].n_candidates; j++) {
                    int s = reads[r].species_indices[j];
                    int L = p->amp_lens[s * M + m_idx];
                    if (L > 0) {
                        gamma_L_sum += gamma[r][j] * (double)L;
                        gamma_sum += gamma[r][j];
                    }
                }
            }
            if (gamma_L_sum > 0) {
                p->lambda = gamma_sum / gamma_L_sum;
                if (p->lambda < 1e-6) p->lambda = 1e-6;
                if (p->lambda > 0.1) p->lambda = 0.1;
            }
        }
    }

    free(eff_counts);
}

/* --- Single EM run --- */
static double em_run_once(em_params_t *p, const em_read_t *reads, int n_reads,
                          const em_config_t *cfg, int single_marker,
                          int *out_iters) {
    int S = p->S;
    (void)S;
    double prev_ll = -INFINITY;
    int iter;

    /* Allocate gamma (sparse: per-read, per-candidate) */
    double **gamma = (double **)hs_malloc((size_t)n_reads * sizeof(double *));
    for (int r = 0; r < n_reads; r++)
        gamma[r] = (double *)hs_calloc((size_t)reads[r].n_candidates, sizeof(double));

    for (iter = 0; iter < cfg->max_iter; iter++) {
        double ll = e_step(p, reads, n_reads, gamma);
        m_step(p, reads, n_reads, gamma, cfg, single_marker);

        double rel_change = fabs(ll - prev_ll) / (fabs(ll) + 1e-10);
        if (iter > 0 && rel_change < cfg->conv_threshold) {
            iter++;
            break;
        }
        prev_ll = ll;
    }

    double final_ll = e_step(p, reads, n_reads, gamma);

    for (int r = 0; r < n_reads; r++) free(gamma[r]);
    free(gamma);

    *out_iters = iter;
    return final_ll;
}

/* --- Public API --- */

em_result_t *em_fit(const em_read_t *reads, int n_reads,
                     int n_species, int n_markers,
                     const int *amplicon_lengths,
                     const em_config_t *config) {
    if (n_reads <= 0 || n_species <= 0) return NULL;

    hs_rng_t rng;
    hs_rng_seed(&rng, config->seed);

    /* Detect single-marker mode: count distinct marker indices across reads */
    int single_marker = 0;
    {
        int marker_seen[8] = {0};  /* HS_MAX_MARKERS = 8 */
        for (int r = 0; r < n_reads; r++) {
            int m = reads[r].marker_idx;
            if (m >= 0 && m < 8) marker_seen[m] = 1;
        }
        int n_markers_seen = 0;
        for (int m = 0; m < 8; m++) n_markers_seen += marker_seen[m];
        if (n_markers_seen <= 1) single_marker = 1;
    }

    if (single_marker) {
        HS_LOG_INFO("Single-marker mode: fixing d from mito CN priors, b=1.0");
    }

    em_params_t *best_p = NULL;
    double best_ll = -INFINITY;
    int best_iters = 0;
    int best_converged = 0;

    int n_restarts = config->n_restarts > 0 ? config->n_restarts : 1;

    for (int restart = 0; restart < n_restarts; restart++) {
        em_params_t *p = params_alloc(n_species, n_markers);
        p->amp_lens = amplicon_lengths;
        p->estimate_degradation = config->estimate_degradation;

        if (restart == 0)
            params_init_uniform(p);
        else
            params_init_random(p, &rng, config);

        /* In single-marker mode, fix d from mito copy numbers and b = 1.0 */
        if (single_marker) {
            if (config->mito_copy_numbers) {
                /* d_s = mito_cn_s / geometric_mean(mito_cn) */
                double log_sum = 0.0;
                int cn_count = 0;
                for (int s = 0; s < n_species; s++) {
                    if (config->mito_copy_numbers[s] > 0) {
                        log_sum += log(config->mito_copy_numbers[s]);
                        cn_count++;
                    }
                }
                double geomean = cn_count > 0 ? exp(log_sum / cn_count) : 1.0;
                for (int s = 0; s < n_species; s++) {
                    p->d[s] = config->mito_copy_numbers[s] > 0
                              ? config->mito_copy_numbers[s] / geomean
                              : 1.0;
                }
            } else {
                for (int s = 0; s < n_species; s++) p->d[s] = 1.0;
            }
            for (int s = 0; s < n_species; s++)
                for (int m = 0; m < n_markers; m++)
                    p->b[s * n_markers + m] = 1.0;
        }

        int iters;
        double ll = em_run_once(p, reads, n_reads, config, single_marker, &iters);

        int converged = (iters < config->max_iter);

        if (ll > best_ll) {
            params_free(best_p);
            best_p = p;
            best_ll = ll;
            best_iters = iters;
            best_converged = converged;
        } else {
            params_free(p);
        }
    }

    /* Build result */
    em_result_t *result = (em_result_t *)hs_calloc(1, sizeof(em_result_t));
    result->n_species = n_species;
    result->n_markers = n_markers;
    result->w = (double *)hs_malloc((size_t)n_species * sizeof(double));
    result->d = (double *)hs_malloc((size_t)n_species * sizeof(double));
    result->b = (double *)hs_malloc((size_t)(n_species * n_markers) * sizeof(double));
    result->w_ci_lo = (double *)hs_calloc((size_t)n_species, sizeof(double));
    result->w_ci_hi = (double *)hs_calloc((size_t)n_species, sizeof(double));

    memcpy(result->w, best_p->w, (size_t)n_species * sizeof(double));
    memcpy(result->d, best_p->d, (size_t)n_species * sizeof(double));
    memcpy(result->b, best_p->b, (size_t)(n_species * n_markers) * sizeof(double));
    result->lambda_proc = best_p->lambda;
    result->log_likelihood = best_ll;
    result->n_iterations = best_iters;
    result->converged = best_converged;

    /* BIC = -2*LL + k*ln(n), k = S-1 (w) + S (d) + S*M (b) [+ 1 (lambda)] */
    int n_params = (n_species - 1) + n_species + n_species * n_markers;
    if (config->estimate_degradation) n_params++;
    result->bic = -2.0 * best_ll + n_params * log((double)n_reads);

    /* Compute confidence intervals */
    if (config->use_advanced_ci)
        em_fisher_info_observed(result, reads, n_reads, n_species, n_markers, amplicon_lengths);
    else
        em_fisher_info(result, reads, n_reads, n_species, n_markers, amplicon_lengths);

    /* Perform Likelihood Ratio Test for species detection */
    if (config->use_full_lrt)
        em_lrt_full(result, reads, n_reads, config, amplicon_lengths);
    else
        em_lrt(result, reads, n_reads, config, amplicon_lengths);

    /* Post-EM species pruning: zero out species below threshold, renormalize */
    if (config->prune_threshold > 0.0) {
        double sum_w = 0.0;
        for (int s = 0; s < n_species; s++) {
            if (result->w[s] < config->prune_threshold)
                result->w[s] = 0.0;
            sum_w += result->w[s];
        }
        if (sum_w > 0.0) {
            for (int s = 0; s < n_species; s++)
                result->w[s] /= sum_w;
        }
    }

    params_free(best_p);
    return result;
}

void em_fisher_info(em_result_t *result,
                     const em_read_t *reads, int n_reads,
                     int n_species, int n_markers,
                     const int *amplicon_lengths) {
    /* Observed Fisher information approximation for w.
     * For Dirichlet-multinomial, Var(w_s) ~ w_s * (1 - w_s) / N_eff
     * where N_eff accounts for classification uncertainty.
     */
    int M = n_markers;

    /* Compute effective sample size per species */
    double *eff_n = (double *)hs_calloc((size_t)n_species, sizeof(double));
    double **gamma = (double **)hs_malloc((size_t)n_reads * sizeof(double *));
    for (int r = 0; r < n_reads; r++) {
        gamma[r] = (double *)hs_calloc((size_t)reads[r].n_candidates, sizeof(double));
        int nc = reads[r].n_candidates;
        int m = reads[r].marker_idx;
        if (m < 0) m = 0;
        double max_log = -INFINITY;
        double *log_g = (double *)hs_malloc((size_t)nc * sizeof(double));
        for (int j = 0; j < nc; j++) {
            int s = reads[r].species_indices[j];
            double c = reads[r].containments[j];
            if (c < 1e-300) c = 1e-300;
            log_g[j] = log(result->w[s]) + log(result->d[s]) +
                       log(result->b[s * M + m]) + log(c);
            if (result->lambda_proc > 1e-8 && amplicon_lengths) {
                int L = amplicon_lengths[s * M + m];
                if (L > 0) log_g[j] -= result->lambda_proc * (double)L;
            }
            if (log_g[j] > max_log) max_log = log_g[j];
        }
        double sum = 0.0;
        for (int j = 0; j < nc; j++) sum += exp(log_g[j] - max_log);
        for (int j = 0; j < nc; j++) {
            gamma[r][j] = exp(log_g[j] - max_log) / sum;
            eff_n[reads[r].species_indices[j]] += gamma[r][j];
        }
        free(log_g);
    }

    /* 95% CI using normal approximation on sqrt(w) transform */
    for (int s = 0; s < n_species; s++) {
        if (eff_n[s] > 0 && result->w[s] > 1e-10) {
            double se = sqrt(result->w[s] * (1.0 - result->w[s]) / (eff_n[s] + 1.0));
            result->w_ci_lo[s] = result->w[s] - 1.96 * se;
            result->w_ci_hi[s] = result->w[s] + 1.96 * se;
            if (result->w_ci_lo[s] < 0.0) result->w_ci_lo[s] = 0.0;
            if (result->w_ci_hi[s] > 1.0) result->w_ci_hi[s] = 1.0;
        } else {
            result->w_ci_lo[s] = 0.0;
            result->w_ci_hi[s] = result->w[s];
        }
    }

    for (int r = 0; r < n_reads; r++) free(gamma[r]);
    free(gamma);
    free(eff_n);
}

/* Simple approximation of Chi-square survival function for df=1 (LRT) */
static double chisq_q_df1(double x) {
    if (x <= 0) return 1.0;
    /* For df=1, Q(x) = erfc(sqrt(x/2)) */
    return erfc(sqrt(x * 0.5));
}

void em_lrt(em_result_t *result,
            const em_read_t *reads, int n_reads,
            const em_config_t *config,
            const int *amplicon_lengths) {
    (void)config;
    int S = result->n_species;
    int M = result->n_markers;
    result->lrt_scores = (double *)hs_calloc((size_t)S, sizeof(double));
    result->p_values = (double *)hs_calloc((size_t)S, sizeof(double));

    /* For each species, compute log-likelihood of null model (w_s = 0) */
    for (int s = 0; s < S; s++) {
        if (result->w[s] < 1e-6) {
            result->lrt_scores[s] = 0.0;
            result->p_values[s] = 1.0;
            continue;
        }

        /* Re-normalize weights without species s */
        double *w_null = (double *)hs_malloc((size_t)S * sizeof(double));
        double sum_null = 0.0;
        for (int i = 0; i < S; i++) {
            w_null[i] = (i == s) ? 0.0 : result->w[i];
            sum_null += w_null[i];
        }
        if (sum_null > 0) {
            for (int i = 0; i < S; i++) w_null[i] /= sum_null;
        } else {
            /* All species nullified? Should not happen if data exists */
            free(w_null);
            continue;
        }

        /* Compute log-likelihood under null model */
        double ll_null = 0.0;
        for (int r = 0; r < n_reads; r++) {
            int nc = reads[r].n_candidates;
            int m = reads[r].marker_idx;
            if (m < 0) m = 0;

            double max_log = -INFINITY;
            double *log_g = (double *)hs_malloc((size_t)nc * sizeof(double));
            int has_valid = 0;

            for (int j = 0; j < nc; j++) {
                int sp = reads[r].species_indices[j];
                if (w_null[sp] <= 0) continue;
                double c = reads[r].containments[j];
                log_g[j] = log(w_null[sp]) + log(result->d[sp]) +
                           log(result->b[sp * M + m]) + log(c);
                if (result->lambda_proc > 1e-8 && amplicon_lengths) {
                    int L = amplicon_lengths[sp * M + m];
                    if (L > 0) log_g[j] -= result->lambda_proc * (double)L;
                }
                if (log_g[j] > max_log) max_log = log_g[j];
                has_valid = 1;
            }

            if (has_valid) {
                double sum = 0.0;
                for (int j = 0; j < nc; j++) {
                    int sp = reads[r].species_indices[j];
                    if (w_null[sp] > 0) sum += exp(log_g[j] - max_log);
                }
                ll_null += max_log + log(sum);
            } else {
                /* Read cannot be explained without species s */
                ll_null += -100.0; /* Penalty for unexplained read */
            }
            free(log_g);
        }

        /* LRT statistic = 2 * (LL_full - LL_null) */
        double lrt = 2.0 * (result->log_likelihood - ll_null);
        if (lrt < 0) lrt = 0; /* numerical precision */
        result->lrt_scores[s] = lrt;
        result->p_values[s] = chisq_q_df1(lrt);

        free(w_null);
    }
}

/* --- Observed Fisher Information CIs (Louis 1982) --- */
void em_fisher_info_observed(em_result_t *result,
                              const em_read_t *reads, int n_reads,
                              int n_species, int n_markers,
                              const int *amplicon_lengths) {
    int S = n_species, M = n_markers;

    /* Re-run E-step with final parameters to get gamma */
    double **gamma = (double **)hs_malloc((size_t)n_reads * sizeof(double *));
    for (int r = 0; r < n_reads; r++) {
        int nc = reads[r].n_candidates;
        gamma[r] = (double *)hs_calloc((size_t)nc, sizeof(double));
        int m = reads[r].marker_idx;
        if (m < 0) m = 0;
        double max_log = -INFINITY;
        double *log_g = (double *)hs_malloc((size_t)nc * sizeof(double));
        for (int j = 0; j < nc; j++) {
            int s = reads[r].species_indices[j];
            double c = reads[r].containments[j];
            if (c < 1e-300) c = 1e-300;
            log_g[j] = log(result->w[s]) + log(result->d[s]) +
                       log(result->b[s * M + m]) + log(c);
            if (result->lambda_proc > 1e-8 && amplicon_lengths) {
                int L = amplicon_lengths[s * M + m];
                if (L > 0) log_g[j] -= result->lambda_proc * (double)L;
            }
            if (log_g[j] > max_log) max_log = log_g[j];
        }
        double sum = 0.0;
        for (int j = 0; j < nc; j++) sum += exp(log_g[j] - max_log);
        for (int j = 0; j < nc; j++)
            gamma[r][j] = exp(log_g[j] - max_log) / sum;
        free(log_g);
    }

    /* Compute observed Fisher information per species with Louis correction */
    for (int s = 0; s < S; s++) {
        if (result->w[s] < 1e-10) {
            result->w_ci_lo[s] = 0.0;
            result->w_ci_hi[s] = 0.0;
            continue;
        }

        double I_complete = 0.0;   /* Complete-data Fisher info */
        double I_missing = 0.0;    /* Missing-data correction */

        for (int r = 0; r < n_reads; r++) {
            /* Find gamma for species s in this read */
            double gamma_rs = 0.0;
            for (int j = 0; j < reads[r].n_candidates; j++) {
                if (reads[r].species_indices[j] == s) {
                    gamma_rs = gamma[r][j];
                    break;
                }
            }

            /* Score contribution: d/dw_s log p(r|theta) */
            double ws = result->w[s];
            double score = gamma_rs / ws - (1.0 - gamma_rs) / (1.0 - ws);

            /* Complete-data Fisher info: sum of score^2 */
            I_complete += score * score;

            /* Louis missing information: Var_Z[score | r, theta] */
            double inv_term = 1.0 / ws + 1.0 / (1.0 - ws);
            I_missing += gamma_rs * (1.0 - gamma_rs) * inv_term * inv_term;
        }

        /* Observed info = complete - missing */
        double I_obs = I_complete - I_missing;
        if (I_obs < 1e-10) I_obs = 1e-10;

        double se = 1.0 / sqrt(I_obs);
        result->w_ci_lo[s] = result->w[s] - 1.96 * se;
        result->w_ci_hi[s] = result->w[s] + 1.96 * se;
        if (result->w_ci_lo[s] < 0.0) result->w_ci_lo[s] = 0.0;
        if (result->w_ci_hi[s] > 1.0) result->w_ci_hi[s] = 1.0;
    }

    for (int r = 0; r < n_reads; r++) free(gamma[r]);
    free(gamma);
}

/* --- Brent's method for lambda optimization --- */
static double brent_obs_ll(double lambda, const em_params_t *p,
                            const em_read_t *reads, int n_reads) {
    /* Compute the observed-data log-likelihood at a given lambda,
     * holding all other parameters (w, d, b) fixed.
     * This directly maximizes the observed LL with respect to lambda. */
    int M = p->M;
    double ll = 0.0;
    for (int r = 0; r < n_reads; r++) {
        int nc = reads[r].n_candidates;
        int m = reads[r].marker_idx;
        if (m < 0) m = 0;

        double max_log = -INFINITY;
        double *log_terms = (double *)hs_malloc((size_t)nc * sizeof(double));
        for (int j = 0; j < nc; j++) {
            int s = reads[r].species_indices[j];
            double c = reads[r].containments[j];
            if (c < 1e-300) c = 1e-300;
            log_terms[j] = log(p->w[s]) + log(p->d[s]) + log(p->b[s * M + m])
                         - lambda * (double)(p->amp_lens ? p->amp_lens[s * M + m] : 0)
                         + log(c);
            if (log_terms[j] > max_log) max_log = log_terms[j];
        }
        double sum_exp = 0.0;
        for (int j = 0; j < nc; j++) sum_exp += exp(log_terms[j] - max_log);
        ll += max_log + log(sum_exp);
        free(log_terms);
    }
    return ll;
}

static double brent_lambda_update(const em_params_t *p, const em_read_t *reads,
                                   int n_reads, double **gamma) {
    (void)gamma; /* Not needed — we maximize observed LL directly */
    /* Brent's method to maximize observed LL w.r.t. lambda over [a, b] */
    double a = 1e-6, b = 0.1;
    double tol = 1e-8;
    int max_iter = 50;

    /* Golden ratio constants */
    double golden = 0.3819660112501051;  /* (3 - sqrt(5)) / 2 */

    double x = a + golden * (b - a);
    double w_br = x, v = x;
    double fx = -brent_obs_ll(x, p, reads, n_reads); /* minimize -Q */
    double fw = fx, fv = fx;
    double d_br = 0.0, e = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        double midpoint = 0.5 * (a + b);
        double tol1 = tol * fabs(x) + 1e-10;
        double tol2 = 2.0 * tol1;

        if (fabs(x - midpoint) <= tol2 - 0.5 * (b - a))
            return x;

        /* Try parabolic interpolation */
        int use_golden = 1;
        double u = 0;
        if (fabs(e) > tol1) {
            double r_br = (x - w_br) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p_br = (x - v) * q - (x - w_br) * r_br;
            q = 2.0 * (q - r_br);
            if (q > 0) p_br = -p_br; else q = -q;

            if (fabs(p_br) < fabs(0.5 * q * e) && p_br > q * (a - x) && p_br < q * (b - x)) {
                d_br = p_br / q;
                u = x + d_br;
                if ((u - a) < tol2 || (b - u) < tol2)
                    d_br = (x < midpoint) ? tol1 : -tol1;
                use_golden = 0;
            }
        }
        if (use_golden) {
            e = (x < midpoint) ? b - x : a - x;
            d_br = golden * e;
        }

        u = (fabs(d_br) >= tol1) ? x + d_br : x + ((d_br > 0) ? tol1 : -tol1);
        double fu = -brent_obs_ll(u, p, reads, n_reads);

        if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w_br; fv = fw;
            w_br = x; fw = fx;
            x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w_br == x) {
                v = w_br; fv = fw;
                w_br = u; fw = fu;
            } else if (fu <= fv || v == x || v == w_br) {
                v = u; fv = fu;
            }
        }
        if (!use_golden) e = d_br;
    }
    return x;
}

/* --- Full nested-model LRT --- */
void em_lrt_full(em_result_t *result,
                 const em_read_t *reads, int n_reads,
                 const em_config_t *config,
                 const int *amplicon_lengths) {
    int S = result->n_species;
    int M = result->n_markers;
    result->lrt_scores = (double *)hs_calloc((size_t)S, sizeof(double));
    result->p_values = (double *)hs_calloc((size_t)S, sizeof(double));

    for (int s_test = 0; s_test < S; s_test++) {
        if (result->w[s_test] < 1e-6) {
            result->lrt_scores[s_test] = 0.0;
            result->p_values[s_test] = 1.0;
            continue;
        }

        /* Build reduced reads: remove species s_test from each read's candidates */
        int n_reduced = 0;
        em_read_t *reduced = (em_read_t *)hs_calloc((size_t)n_reads, sizeof(em_read_t));

        /* Species index remapping: shift indices > s_test down by 1 */
        for (int r = 0; r < n_reads; r++) {
            int nc_orig = reads[r].n_candidates;
            int nc_new = 0;

            /* Count candidates excluding s_test */
            for (int j = 0; j < nc_orig; j++)
                if (reads[r].species_indices[j] != s_test) nc_new++;

            if (nc_new == 0) continue; /* Read only had s_test — drop it */

            reduced[n_reduced].marker_idx = reads[r].marker_idx;
            reduced[n_reduced].n_candidates = nc_new;
            reduced[n_reduced].species_indices = (int *)hs_malloc((size_t)nc_new * sizeof(int));
            reduced[n_reduced].containments = (double *)hs_malloc((size_t)nc_new * sizeof(double));

            int k = 0;
            for (int j = 0; j < nc_orig; j++) {
                int sp = reads[r].species_indices[j];
                if (sp == s_test) continue;
                /* Remap: indices > s_test shift down by 1 */
                reduced[n_reduced].species_indices[k] = sp > s_test ? sp - 1 : sp;
                reduced[n_reduced].containments[k] = reads[r].containments[j];
                k++;
            }
            n_reduced++;
        }

        /* Build reduced amplicon lengths (remove row s_test) */
        int S_red = S - 1;
        int *amp_red = NULL;
        if (amplicon_lengths) {
            amp_red = (int *)hs_calloc((size_t)(S_red * M), sizeof(int));
            int dst = 0;
            for (int s = 0; s < S; s++) {
                if (s == s_test) continue;
                for (int m = 0; m < M; m++)
                    amp_red[dst * M + m] = amplicon_lengths[s * M + m];
                dst++;
            }
        }

        /* Re-fit EM with S-1 species */
        em_config_t reduced_cfg = *config;
        reduced_cfg.n_restarts = 1;        /* Warm-start sufficient */
        reduced_cfg.max_iter = 50;         /* Fewer iterations needed */
        reduced_cfg.use_advanced_ci = 0;   /* Disable advanced methods in sub-fit */
        reduced_cfg.use_brent_lambda = 0;
        reduced_cfg.use_full_lrt = 0;
        reduced_cfg.prune_threshold = 0.0;
        reduced_cfg.mito_copy_numbers = NULL; /* Not needed for sub-fit */

        double ll_reduced = -INFINITY;
        if (n_reduced > 0 && S_red > 0) {
            em_result_t *res_red = em_fit(reduced, n_reduced, S_red, M,
                                           amp_red, &reduced_cfg);
            if (res_red) {
                ll_reduced = res_red->log_likelihood;
                em_result_destroy(res_red);
            }
        }

        /* LRT = 2 * (LL_full - LL_reduced) */
        double lrt = 2.0 * (result->log_likelihood - ll_reduced);
        if (lrt < 0) lrt = 0;
        result->lrt_scores[s_test] = lrt;
        result->p_values[s_test] = chisq_q_df1(lrt);

        /* Cleanup */
        free(amp_red);
        for (int r = 0; r < n_reduced; r++) {
            free(reduced[r].species_indices);
            free(reduced[r].containments);
        }
        free(reduced);
    }
}

em_read_t *em_reads_from_classify(const void *results_ptr, int n_reads,
                                   int *out_n_em_reads) {
    const read_result_t *results = (const read_result_t *)results_ptr;
    /* Count classified reads */
    int n_em = 0;
    for (int i = 0; i < n_reads; i++)
        if (results[i].is_classified && results[i].n_hits > 0) n_em++;

    em_read_t *em_reads = (em_read_t *)hs_calloc((size_t)n_em, sizeof(em_read_t));
    int idx = 0;
    for (int i = 0; i < n_reads; i++) {
        if (!results[i].is_classified || results[i].n_hits == 0) continue;
        em_reads[idx].marker_idx = results[i].marker_idx;
        em_reads[idx].n_candidates = results[i].n_hits;
        em_reads[idx].species_indices = (int *)hs_malloc(
            (size_t)results[i].n_hits * sizeof(int));
        em_reads[idx].containments = (double *)hs_malloc(
            (size_t)results[i].n_hits * sizeof(double));
        for (int j = 0; j < results[i].n_hits; j++) {
            em_reads[idx].species_indices[j] = results[i].hits[j].species_idx;
            em_reads[idx].containments[j] = results[i].hits[j].containment;
        }
        idx++;
    }
    *out_n_em_reads = n_em;
    return em_reads;
}

void em_reads_free(em_read_t *reads, int n) {
    if (!reads) return;
    for (int i = 0; i < n; i++) {
        free(reads[i].species_indices);
        free(reads[i].containments);
    }
    free(reads);
}

void em_result_destroy(em_result_t *r) {
    if (!r) return;
    free(r->w); free(r->d); free(r->b);
    free(r->w_ci_lo); free(r->w_ci_hi);
    free(r->lrt_scores); free(r->p_values);
    free(r);
}
