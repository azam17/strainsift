#ifndef SPECIESID_EM_H
#define SPECIESID_EM_H

#include <stdint.h>

typedef struct {
    int max_iter;
    double conv_threshold;
    int n_restarts;
    double alpha;              /* Dirichlet prior on w */
    double d_mu, d_sigma;      /* LogNormal prior on DNA yield */
    double b_mu, b_sigma;      /* LogNormal prior on PCR bias */
    int estimate_degradation;
    double min_containment;
    uint64_t seed;
    double *mito_copy_numbers; /* [S] mito CN per species (NULL = estimate d freely) */
    double prune_threshold;    /* Post-EM pruning: remove species below this weight (0 = disabled) */
    int use_advanced_ci;       /* 0 = Wald (default), 1 = observed Fisher information */
    int use_brent_lambda;      /* 0 = closed-form (default), 1 = Brent's method */
    int use_full_lrt;          /* 0 = profile LRT (default), 1 = full nested-model refit */
} em_config_t;

typedef struct {
    double *w;                 /* [S] species weight fractions */
    double *d;                 /* [S] DNA yield factors */
    double *b;                 /* [S * M] PCR bias matrix (row-major) */
    double lambda_proc;        /* Degradation rate */
    double *w_ci_lo, *w_ci_hi; /* [S] 95% CIs on weight fractions */
    double log_likelihood;
    double bic;
    int n_iterations;
    int converged;
    double *lrt_scores;   /* [S] Likelihood ratio test scores */
    double *p_values;     /* [S] LRT p-values */
    int n_species;
    int n_markers;
} em_result_t;

/* Sparse read data for EM */
typedef struct {
    int marker_idx;
    int *species_indices;      /* Candidate species for this read */
    double *containments;      /* Containment scores */
    int n_candidates;
} em_read_t;

em_config_t em_config_default(void);

/* THE CORE: fit multi-marker bias-corrected EM */
em_result_t *em_fit(const em_read_t *reads, int n_reads,
                     int n_species, int n_markers,
                     const int *amplicon_lengths,
                     const em_config_t *config);

/* Compute Fisher information for confidence intervals (Wald approximation) */
void em_fisher_info(em_result_t *result,
                     const em_read_t *reads, int n_reads,
                     int n_species, int n_markers,
                     const int *amplicon_lengths);

/* Observed Fisher information CIs with Louis (1982) missing-data correction */
void em_fisher_info_observed(em_result_t *result,
                              const em_read_t *reads, int n_reads,
                              int n_species, int n_markers,
                              const int *amplicon_lengths);

/* Perform Likelihood Ratio Test (LRT) for species presence (profile) */
void em_lrt(em_result_t *result,
            const em_read_t *reads, int n_reads,
            const em_config_t *config,
            const int *amplicon_lengths);

/* Full nested-model LRT: re-fits EM with each species removed */
void em_lrt_full(em_result_t *result,
                 const em_read_t *reads, int n_reads,
                 const em_config_t *config,
                 const int *amplicon_lengths);

/* Build em_read_t array from classify results */
em_read_t *em_reads_from_classify(const void *results, int n_reads,
                                   int *out_n_em_reads);

void em_reads_free(em_read_t *reads, int n);
void em_result_destroy(em_result_t *r);

#endif /* SPECIESID_EM_H */
