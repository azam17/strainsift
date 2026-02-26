#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "em.h"
#include "calibrate.h"
#include "utils.h"

static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "  FAIL: %s (line %d)\n", msg, __LINE__); \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

#define ASSERT_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        fprintf(stderr, "  FAIL: %s: %.6f != %.6f (tol=%.6f) (line %d)\n", \
                msg, (a), (b), (tol), __LINE__); \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

/* Helper: generate synthetic EM reads for a 2-species mixture */
static em_read_t *make_reads_2species(int n_reads, int n_markers,
                                       double w0, double w1,
                                       double *bias, /* [2*M] */
                                       uint64_t seed, int *out_n) {
    hs_rng_t rng;
    hs_rng_seed(&rng, seed);

    em_read_t *reads = (em_read_t *)hs_calloc((size_t)n_reads, sizeof(em_read_t));
    int count = 0;

    for (int r = 0; r < n_reads; r++) {
        int m = r % n_markers;
        /* Compute probabilities */
        double p0 = w0 * (bias ? bias[0 * n_markers + m] : 1.0);
        double p1 = w1 * (bias ? bias[1 * n_markers + m] : 1.0);
        double total = p0 + p1;
        p0 /= total; p1 /= total;

        /* Assign read */
        double u = hs_rng_uniform(&rng);
        int true_sp = (u < p0) ? 0 : 1;

        reads[count].marker_idx = m;
        reads[count].n_candidates = 2;
        reads[count].species_indices = (int *)hs_malloc(2 * sizeof(int));
        reads[count].containments = (double *)hs_malloc(2 * sizeof(double));
        reads[count].species_indices[0] = 0;
        reads[count].species_indices[1] = 1;
        /* Simulate containment: high for true species, low for false */
        if (true_sp == 0) {
            reads[count].containments[0] = 0.8 + hs_rng_uniform(&rng) * 0.2;
            reads[count].containments[1] = hs_rng_uniform(&rng) * 0.2;
        } else {
            reads[count].containments[0] = hs_rng_uniform(&rng) * 0.2;
            reads[count].containments[1] = 0.8 + hs_rng_uniform(&rng) * 0.2;
        }
        count++;
    }
    *out_n = count;
    return reads;
}

static void test_em_basic_50_50(void) {
    printf("  test_em_basic_50_50...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.5, 0.5, NULL, 42, &n_reads);

    em_config_t cfg = em_config_default();
    cfg.n_restarts = 3;
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");
    ASSERT(res->converged, "Marked as converged");

    /* Both weights should be near 0.5 */
    ASSERT_NEAR(res->w[0], 0.5, 0.1, "w[0] near 0.5");
    ASSERT_NEAR(res->w[1], 0.5, 0.1, "w[1] near 0.5");
    ASSERT(res->w[0] + res->w[1] > 0.99, "Weights sum to ~1");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_90_10(void) {
    printf("  test_em_90_10...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.9, 0.1, NULL, 123, &n_reads);

    em_config_t cfg = em_config_default();
    cfg.n_restarts = 3;
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");

    /* Dominant species should have w > 0.7 */
    double dominant = res->w[0] > res->w[1] ? res->w[0] : res->w[1];
    double minor = res->w[0] < res->w[1] ? res->w[0] : res->w[1];
    ASSERT(dominant > 0.7, "Dominant species w > 0.7");
    ASSERT(minor < 0.3, "Minor species w < 0.3");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_bias_recovery(void) {
    printf("  test_em_bias_recovery...\n");
    /* 50/50 mixture with 3x PCR bias on species 1 at marker 0 */
    double bias[6] = { 1.0, 1.0, 1.0,   /* species 0: no bias */
                       3.0, 1.0, 1.0 };  /* species 1: 3x at COI */

    int n_reads;
    em_read_t *reads = make_reads_2species(6000, 3, 0.5, 0.5, bias, 456, &n_reads);

    em_config_t cfg = em_config_default();
    cfg.n_restarts = 5;
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged with bias");

    /* Despite biased reads, EM should recover near 50/50 weights */
    ASSERT_NEAR(res->w[0], 0.5, 0.15, "w[0] near 0.5 despite bias");
    ASSERT_NEAR(res->w[1], 0.5, 0.15, "w[1] near 0.5 despite bias");

    /* Bias should be recovered in b matrix */
    /* Species 1, marker 0 should have higher bias than marker 1 */
    ASSERT(res->b[1 * 3 + 0] > res->b[1 * 3 + 1] * 0.5,
           "Bias pattern partially recovered");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_confidence_intervals(void) {
    printf("  test_em_confidence_intervals...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.5, 0.5, NULL, 789, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");

    /* CIs should contain the true value (0.5) */
    ASSERT(res->w_ci_lo[0] < 0.5 && res->w_ci_hi[0] > 0.5,
           "CI[0] contains true value");
    ASSERT(res->w_ci_lo[1] < 0.5 && res->w_ci_hi[1] > 0.5,
           "CI[1] contains true value");

    /* CIs should be reasonable width */
    double width0 = res->w_ci_hi[0] - res->w_ci_lo[0];
    ASSERT(width0 < 0.5, "CI width < 0.5");
    ASSERT(width0 > 0.001, "CI width > 0.001");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_bic(void) {
    printf("  test_em_bic...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(1000, 3, 0.5, 0.5, NULL, 321, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");
    ASSERT(isfinite(res->bic), "BIC is finite");
    ASSERT(res->log_likelihood < 0, "Log-likelihood is negative");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_single_species(void) {
    printf("  test_em_single_species...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(1000, 3, 1.0, 0.0, NULL, 654, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");
    /* Dominant species should get almost all weight */
    double dominant = res->w[0] > res->w[1] ? res->w[0] : res->w[1];
    ASSERT(dominant > 0.8, "Single species gets most weight");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_d_factors(void) {
    printf("  test_em_d_factors...\n");
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.5, 0.5, NULL, 111, &n_reads);

    em_config_t cfg = em_config_default();
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    em_result_t *res = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg);
    ASSERT(res != NULL, "EM converged");

    /* DNA yield factors should be near 1.0 for equal species */
    ASSERT(res->d[0] > 0.3 && res->d[0] < 3.0, "d[0] reasonable range");
    ASSERT(res->d[1] > 0.3 && res->d[1] < 3.0, "d[1] reasonable range");

    em_result_destroy(res);
    em_reads_free(reads, n_reads);
}

static void test_em_single_marker_mode(void) {
    printf("  test_em_single_marker_mode...\n");
    /* Simulate a single-marker (16S) scenario:
     * Species 0 has 2x mito CN vs species 1, so DNA yield is biased.
     * True composition is 50/50, but without mito CN correction the EM
     * would over-estimate the high-CN species from raw read counts. */
    int n_total = 3000;
    int n_reads = 0;
    em_read_t *reads = (em_read_t *)hs_calloc((size_t)n_total, sizeof(em_read_t));

    hs_rng_t rng;
    hs_rng_seed(&rng, 999);

    /* Species 0: mito_cn=2000, species 1: mito_cn=1000
     * With 50/50 true weight, species 0 produces 2x more reads */
    double mito_cn[2] = { 2000.0, 1000.0 };
    double eff0 = 0.5 * mito_cn[0];
    double eff1 = 0.5 * mito_cn[1];
    double total_eff = eff0 + eff1;
    double p0 = eff0 / total_eff;  /* ~0.667 of reads from sp0 */

    for (int r = 0; r < n_total; r++) {
        double u = hs_rng_uniform(&rng);
        int true_sp = (u < p0) ? 0 : 1;

        reads[n_reads].marker_idx = 0;  /* single marker */
        reads[n_reads].n_candidates = 2;
        reads[n_reads].species_indices = (int *)hs_malloc(2 * sizeof(int));
        reads[n_reads].containments = (double *)hs_malloc(2 * sizeof(double));
        reads[n_reads].species_indices[0] = 0;
        reads[n_reads].species_indices[1] = 1;
        if (true_sp == 0) {
            reads[n_reads].containments[0] = 0.85 + hs_rng_uniform(&rng) * 0.15;
            reads[n_reads].containments[1] = hs_rng_uniform(&rng) * 0.15;
        } else {
            reads[n_reads].containments[0] = hs_rng_uniform(&rng) * 0.15;
            reads[n_reads].containments[1] = 0.85 + hs_rng_uniform(&rng) * 0.15;
        }
        n_reads++;
    }

    /* Without mito CN correction: d is free, should converge near raw proportions */
    em_config_t cfg_no_cn = em_config_default();
    cfg_no_cn.mito_copy_numbers = NULL;
    int amp_lens[2] = { 460, 460 };

    em_result_t *res_no = em_fit(reads, n_reads, 2, 1, amp_lens, &cfg_no_cn);
    ASSERT(res_no != NULL, "EM (no CN) converged");
    /* Without correction, species 0 gets inflated weight (~0.6+) */
    double w0_no = res_no->w[0] > res_no->w[1] ? res_no->w[0] : res_no->w[1];
    ASSERT(w0_no > 0.55, "Without mito CN: dominant species overestimated");

    /* With mito CN correction: d fixed from known CN values */
    em_config_t cfg_cn = em_config_default();
    cfg_cn.mito_copy_numbers = mito_cn;

    em_result_t *res_cn = em_fit(reads, n_reads, 2, 1, amp_lens, &cfg_cn);
    ASSERT(res_cn != NULL, "EM (with CN) converged");
    /* With correction, both weights should be closer to 0.5 */
    ASSERT_NEAR(res_cn->w[0], 0.5, 0.15, "With mito CN: w[0] near 0.5");
    ASSERT_NEAR(res_cn->w[1], 0.5, 0.15, "With mito CN: w[1] near 0.5");

    /* d should reflect mito CN ratio (2000/1000 -> d[0]/d[1] ~ 2) */
    double d_ratio = res_cn->d[0] / res_cn->d[1];
    ASSERT(d_ratio > 1.5 && d_ratio < 2.5,
           "d ratio reflects mito CN ratio (~2)");

    /* b should be 1.0 in single-marker mode */
    ASSERT_NEAR(res_cn->b[0], 1.0, 0.01, "b[0]=1.0 in single-marker mode");
    ASSERT_NEAR(res_cn->b[1], 1.0, 0.01, "b[1]=1.0 in single-marker mode");

    em_result_destroy(res_no);
    em_result_destroy(res_cn);
    em_reads_free(reads, n_reads);
}

static void test_em_degradation(void) {
    printf("  test_em_degradation...\n");

    /* Two species, two markers with different amplicon lengths.
     * Species 0: marker 0 = 350bp, marker 1 = 700bp
     * Species 1: marker 0 = 350bp, marker 1 = 700bp
     * With degradation, shorter amplicons should get more reads.
     */
    int n_species = 2, n_markers = 2;
    int amp_lens[4] = {350, 700, 350, 700}; /* [S*M]: s0m0, s0m1, s1m0, s1m1 */
    double true_lambda = 0.002;

    int n_reads = 800;
    em_read_t *reads = (em_read_t *)hs_calloc((size_t)n_reads, sizeof(em_read_t));
    hs_rng_t rng;
    hs_rng_seed(&rng, 123);

    for (int i = 0; i < n_reads; i++) {
        int true_sp = (hs_rng_uniform(&rng) < 0.7) ? 0 : 1;
        double surv0 = exp(-true_lambda * 350.0);
        double surv1 = exp(-true_lambda * 700.0);
        double p0 = surv0 / (surv0 + surv1);
        int marker = (hs_rng_uniform(&rng) < p0) ? 0 : 1;

        reads[i].marker_idx = marker;
        reads[i].n_candidates = n_species;
        reads[i].species_indices = (int *)hs_malloc((size_t)n_species * sizeof(int));
        reads[i].containments = (double *)hs_malloc((size_t)n_species * sizeof(double));

        for (int s = 0; s < n_species; s++) {
            reads[i].species_indices[s] = s;
            reads[i].containments[s] = (s == true_sp) ? 0.95 : 0.05;
        }
    }

    /* Run WITH degradation correction */
    em_config_t cfg = em_config_default();
    cfg.estimate_degradation = 1;

    em_result_t *r = em_fit(reads, n_reads, n_species, n_markers, amp_lens, &cfg);
    ASSERT(r != NULL, "EM result not NULL");
    ASSERT(r->lambda_proc > 1e-6, "Lambda estimated > 0");
    ASSERT_NEAR(r->w[0], 0.7, 0.15, "Degradation-corrected weight near 0.7");

    /* BIC should include lambda parameter */
    int n_params_expected = (n_species - 1) + n_species + n_species * n_markers + 1;
    double bic_expected = -2.0 * r->log_likelihood + n_params_expected * log((double)n_reads);
    ASSERT_NEAR(r->bic, bic_expected, 0.1, "BIC includes lambda parameter");

    em_result_destroy(r);
    em_reads_free(reads, n_reads);
}

static void test_em_no_degradation_regression(void) {
    printf("  test_em_no_degradation_regression...\n");

    int n_reads;
    em_read_t *reads = make_reads_2species(500, 1, 0.5, 0.5, NULL, 777, &n_reads);
    int amp_lens[2] = {400, 400};

    em_config_t cfg = em_config_default();
    cfg.estimate_degradation = 0;

    em_result_t *r = em_fit(reads, n_reads, 2, 1, amp_lens, &cfg);
    ASSERT(r != NULL, "EM result not NULL");
    ASSERT_NEAR(r->w[0], 0.5, 0.15, "50/50 mixture recovered");

    /* BIC should NOT include lambda */
    int n_params_expected = (2 - 1) + 2 + 2 * 1;
    double bic_expected = -2.0 * r->log_likelihood + n_params_expected * log((double)n_reads);
    ASSERT_NEAR(r->bic, bic_expected, 0.1, "BIC without lambda");

    em_result_destroy(r);
    em_reads_free(reads, n_reads);
}

static void test_em_pruning(void) {
    printf("  test_em_pruning...\n");
    /* Pure species 0 input: after pruning, species 0 should be ~100% */
    int n_reads;
    em_read_t *reads = make_reads_2species(1000, 1, 1.0, 0.0, NULL, 888, &n_reads);

    em_config_t cfg = em_config_default();
    cfg.prune_threshold = 0.05;  /* Prune species below 5% */
    int amp_lens[2] = { 460, 460 };

    em_result_t *res = em_fit(reads, n_reads, 2, 1, amp_lens, &cfg);
    ASSERT(res != NULL, "EM with pruning converged");

    /* After pruning, the minor species should be exactly 0 */
    double minor = res->w[0] < res->w[1] ? res->w[0] : res->w[1];
    ASSERT_NEAR(minor, 0.0, 1e-10, "Pruned minor species is 0");

    /* Dominant species should be exactly 1.0 after renormalization */
    double dominant = res->w[0] > res->w[1] ? res->w[0] : res->w[1];
    ASSERT_NEAR(dominant, 1.0, 1e-10, "Dominant species is 1.0 after pruning");

    /* Without pruning, minor species would have some small weight */
    em_config_t cfg_no_prune = em_config_default();
    em_result_t *res2 = em_fit(reads, n_reads, 2, 1, amp_lens, &cfg_no_prune);
    ASSERT(res2 != NULL, "EM without pruning converged");
    double minor2 = res2->w[0] < res2->w[1] ? res2->w[0] : res2->w[1];
    ASSERT(minor2 > 0.0, "Without pruning, minor species has nonzero weight");

    em_result_destroy(res);
    em_result_destroy(res2);
    em_reads_free(reads, n_reads);
}

static void test_calibration_roundtrip(void) {
    printf("  test_calibration_roundtrip...\n");

    calibration_result_t orig;
    orig.d_mu = 0.123456;
    orig.d_sigma = 0.456789;
    orig.b_mu = -0.234567;
    orig.b_sigma = 0.345678;
    orig.n_samples = 5;

    const char *path = "/tmp/test_calibration.cal";
    int rc = calibrate_save(&orig, path);
    ASSERT(rc == 0, "Calibration save succeeded");

    calibration_result_t *loaded = calibrate_load(path);
    ASSERT(loaded != NULL, "Calibration load succeeded");
    ASSERT_NEAR(loaded->d_mu, orig.d_mu, 1e-6, "d_mu roundtrip");
    ASSERT_NEAR(loaded->d_sigma, orig.d_sigma, 1e-6, "d_sigma roundtrip");
    ASSERT_NEAR(loaded->b_mu, orig.b_mu, 1e-6, "b_mu roundtrip");
    ASSERT_NEAR(loaded->b_sigma, orig.b_sigma, 1e-6, "b_sigma roundtrip");
    ASSERT(loaded->n_samples == orig.n_samples, "n_samples roundtrip");

    calibrate_result_destroy(loaded);
    remove(path);
}

static void test_calibration_estimate(void) {
    printf("  test_calibration_estimate...\n");

    calibration_sample_t sample;
    sample.n_species = 2;
    sample.n_markers = 2;
    sample.true_w = (double *)calloc(2, sizeof(double));
    sample.obs_reads = (double *)calloc(4, sizeof(double));
    sample.true_w[0] = 0.7;
    sample.true_w[1] = 0.3;
    sample.obs_reads[0] = 700;  sample.obs_reads[1] = 600;
    sample.obs_reads[2] = 250;  sample.obs_reads[3] = 350;

    calibration_result_t *r = calibrate_estimate(&sample, 1);
    ASSERT(r != NULL, "Calibration estimate not NULL");
    ASSERT(r->d_sigma >= 0.1, "d_sigma >= 0.1");
    ASSERT(r->b_sigma >= 0.1, "b_sigma >= 0.1");

    calibrate_result_destroy(r);
    free(sample.true_w);
    free(sample.obs_reads);
}

/* --- Advanced inference comparison tests --- */

static void test_em_fisher_vs_wald(void) {
    printf("  test_em_fisher_vs_wald...\n");
    /* 90:10 mixture — near-boundary case where Fisher should differ from Wald */
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.9, 0.1, NULL, 555, &n_reads);
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    /* Standard (Wald) */
    em_config_t cfg_wald = em_config_default();
    cfg_wald.n_restarts = 3;
    em_result_t *res_wald = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg_wald);

    /* Advanced (Fisher) */
    em_config_t cfg_fisher = em_config_default();
    cfg_fisher.n_restarts = 3;
    cfg_fisher.use_advanced_ci = 1;
    em_result_t *res_fisher = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg_fisher);

    ASSERT(res_wald != NULL && res_fisher != NULL, "Both methods converged");

    /* Weights should be nearly identical (same EM, just different CIs) */
    ASSERT_NEAR(res_wald->w[0], res_fisher->w[0], 0.01,
                "Weights match between Wald and Fisher");

    /* Fisher CIs should be valid (contain true value 0.9 or 0.1) */
    double dominant_fisher = res_fisher->w[0] > res_fisher->w[1] ? 0 : 1;
    int dom_idx = (int)dominant_fisher;
    /* The true dominant value is 0.9 */
    double true_dom = 0.9;
    double true_min = 0.1;
    if (res_fisher->w[0] < res_fisher->w[1]) {
        dom_idx = 1;
    } else {
        dom_idx = 0;
    }
    int min_idx = 1 - dom_idx;

    ASSERT(res_fisher->w_ci_lo[dom_idx] < true_dom &&
           res_fisher->w_ci_hi[dom_idx] > true_dom,
           "Fisher CI covers true dominant value");
    ASSERT(res_fisher->w_ci_lo[min_idx] < true_min &&
           res_fisher->w_ci_hi[min_idx] > true_min,
           "Fisher CI covers true minor value");

    /* Fisher CIs should differ from Wald CIs (not necessarily tighter in all cases,
     * but they should be different due to the Louis correction) */
    double wald_width_0 = res_wald->w_ci_hi[0] - res_wald->w_ci_lo[0];
    double fisher_width_0 = res_fisher->w_ci_hi[0] - res_fisher->w_ci_lo[0];
    ASSERT(fabs(wald_width_0 - fisher_width_0) > 1e-6 || wald_width_0 < 0.001,
           "Fisher and Wald CIs differ (or both very tight)");

    em_result_destroy(res_wald);
    em_result_destroy(res_fisher);
    em_reads_free(reads, n_reads);
}

static void test_em_brent_vs_closedform(void) {
    printf("  test_em_brent_vs_closedform...\n");
    /* With degradation enabled, compare LL between closed-form and Brent's lambda */
    int n_species = 2, n_markers = 2;
    int amp_lens[4] = {350, 700, 350, 700};
    int n_reads = 1000;

    em_read_t *reads = (em_read_t *)hs_calloc((size_t)n_reads, sizeof(em_read_t));
    hs_rng_t rng;
    hs_rng_seed(&rng, 777);

    double true_lambda = 0.002;
    for (int i = 0; i < n_reads; i++) {
        int true_sp = (hs_rng_uniform(&rng) < 0.7) ? 0 : 1;
        double surv0 = exp(-true_lambda * 350.0);
        double surv1 = exp(-true_lambda * 700.0);
        double p0 = surv0 / (surv0 + surv1);
        int marker = (hs_rng_uniform(&rng) < p0) ? 0 : 1;

        reads[i].marker_idx = marker;
        reads[i].n_candidates = n_species;
        reads[i].species_indices = (int *)hs_malloc((size_t)n_species * sizeof(int));
        reads[i].containments = (double *)hs_malloc((size_t)n_species * sizeof(double));
        for (int s = 0; s < n_species; s++) {
            reads[i].species_indices[s] = s;
            reads[i].containments[s] = (s == true_sp) ? 0.95 : 0.05;
        }
    }

    /* Closed-form lambda */
    em_config_t cfg_cf = em_config_default();
    cfg_cf.estimate_degradation = 1;
    cfg_cf.n_restarts = 3;
    em_result_t *res_cf = em_fit(reads, n_reads, n_species, n_markers, amp_lens, &cfg_cf);

    /* Brent's lambda */
    em_config_t cfg_br = em_config_default();
    cfg_br.estimate_degradation = 1;
    cfg_br.use_brent_lambda = 1;
    cfg_br.n_restarts = 3;
    em_result_t *res_br = em_fit(reads, n_reads, n_species, n_markers, amp_lens, &cfg_br);

    ASSERT(res_cf != NULL && res_br != NULL, "Both methods converged");

    /* Brent should achieve >= log-likelihood of closed-form */
    ASSERT(res_br->log_likelihood >= res_cf->log_likelihood - 1.0,
           "Brent LL >= closed-form LL (within tolerance)");

    /* Both should recover reasonable weights */
    ASSERT_NEAR(res_cf->w[0], 0.7, 0.15, "Closed-form weight near 0.7");
    ASSERT_NEAR(res_br->w[0], 0.7, 0.15, "Brent weight near 0.7");

    em_result_destroy(res_cf);
    em_result_destroy(res_br);
    em_reads_free(reads, n_reads);
}

static void test_em_full_lrt_vs_profile(void) {
    printf("  test_em_full_lrt_vs_profile...\n");
    /* 90:10 mixture — species 0 (dominant) and species 1 (minor) both present.
     * Both LRT methods should detect both as present (p < 0.05).
     * Full LRT should generally have at least as much power. */
    int n_reads;
    em_read_t *reads = make_reads_2species(3000, 3, 0.9, 0.1, NULL, 333, &n_reads);
    int amp_lens[6] = { 658, 425, 560, 658, 425, 560 };

    /* Profile LRT (standard) */
    em_config_t cfg_profile = em_config_default();
    cfg_profile.n_restarts = 3;
    em_result_t *res_prof = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg_profile);

    /* Full nested-model LRT */
    em_config_t cfg_full = em_config_default();
    cfg_full.n_restarts = 3;
    cfg_full.use_full_lrt = 1;
    em_result_t *res_full = em_fit(reads, n_reads, 2, 3, amp_lens, &cfg_full);

    ASSERT(res_prof != NULL && res_full != NULL, "Both LRT methods converged");

    /* Both should detect both species as present (p < 0.05) */
    ASSERT(res_prof->p_values[0] < 0.05, "Profile LRT: dominant species significant");
    ASSERT(res_prof->p_values[1] < 0.05, "Profile LRT: minor species significant");
    ASSERT(res_full->p_values[0] < 0.05, "Full LRT: dominant species significant");
    ASSERT(res_full->p_values[1] < 0.05, "Full LRT: minor species significant");

    /* Full LRT scores should be non-negative */
    ASSERT(res_full->lrt_scores[0] >= 0, "Full LRT score[0] >= 0");
    ASSERT(res_full->lrt_scores[1] >= 0, "Full LRT score[1] >= 0");

    /* Full LRT p-values should be <= profile p-values (more power from proper refit)
     * Allow small numerical tolerance */
    ASSERT(res_full->p_values[0] <= res_prof->p_values[0] + 0.01,
           "Full LRT p[0] <= profile p[0] (more power)");

    em_result_destroy(res_prof);
    em_result_destroy(res_full);
    em_reads_free(reads, n_reads);
}

int main(void) {
    printf("=== test_em ===\n");
    test_em_basic_50_50();
    test_em_90_10();
    test_em_bias_recovery();
    test_em_confidence_intervals();
    test_em_bic();
    test_em_single_species();
    test_em_d_factors();
    test_em_single_marker_mode();
    test_em_degradation();
    test_em_no_degradation_regression();
    test_em_pruning();
    test_calibration_roundtrip();
    test_calibration_estimate();
    test_em_fisher_vs_wald();
    test_em_brent_vs_closedform();
    test_em_full_lrt_vs_profile();
    printf("=== %d passed, %d failed ===\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
