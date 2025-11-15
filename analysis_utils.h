/* analysis_utils.h - TRTS Analysis Utilities
 *
 * Provides in-memory analysis of TRTS runs using observer callbacks.
 * Blueprint reference: analysis_utils.h/analysis_utils.c from Section 1.
 *
 * Uses simulate_stream with an in-memory observer to collect statistics
 * without file I/O. All metrics are computed from raw rational values
 * using Welford's algorithm for numerical stability.
 *
 * CRITICAL: Analysis functions never modify propagation or write results
 * back into the engine. All evaluation is read-only.
 */

#ifndef TRTS_ANALYSIS_UTILS_H
#define TRTS_ANALYSIS_UTILS_H

#include "config.h"
#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>

/* Run summary containing all collected statistics */
typedef struct {
    /* Final ratio */
    // FIX: Changed mpq_t to Rational for consistency with rational.h functions
    Rational final_ratio;                  /* Last observed υ/β ratio */
    bool ratio_defined;                 /* True if ratio was computable */
    char final_ratio_str[128];          /* String representation "N/D" */
    
    /* Convergence analysis */
    char closest_constant[32];          /* Name of nearest known constant */
    double closest_delta;               /* Distance to nearest constant */
    size_t convergence_tick;            /* First tick within 1e-5 of constant */
    
    /* Classification */
    char pattern[32];                   /* "divergent", "fixed point", etc */
    char classification[64];            /* Detailed classification string */
    
    /* Stack statistics */
    char stack_summary[128];            /* Summary of κ stack contents */
    size_t stack_max_depth;             /* Maximum κ stack depth reached */
    
    /* Engine statistics */
    size_t engine_step_count;           /* Number of successful engine steps */
    size_t koppa_sample_count;          /* Number of koppa samples taken */
    
    /* Psi firing statistics */
    size_t psi_fire_count;              /* Total number of ψ fires */
    size_t psi_triple_count;            /* Number of triple ψ fires */
    double psi_spacing_mean;            /* Mean microticks between ψ fires */
    double psi_spacing_stddev;          /* Std dev of spacing */
    
    /* Ratio statistics */
    double ratio_variance;              /* Variance of υ/β samples */
    double ratio_range;                 /* Max - min of samples */
    double ratio_mean;                  /* Mean of υ/β samples */
    double ratio_stddev;                /* Std dev of υ/β samples */
} RunSummary;

/* Initialize run summary (allocate GMP resources) */
void run_summary_init(RunSummary *summary);

/* Clear run summary (free GMP resources) */
void run_summary_clear(RunSummary *summary);

/* Copy run summary */
void run_summary_copy(RunSummary *dest, const RunSummary *src);

/* Analyze latest run using in-memory observer
 *
 * Runs simulation with observer callback, collects all statistics,
 * and populates RunSummary. Uses Welford's algorithm for stable
 * mean/variance computation.
 *
 * Returns: true on success, false on error
 */
bool analyze_latest_run(const Config *config, RunSummary *summary);

/* Convenience wrapper: simulate and analyze in one call */
bool simulate_and_analyze(const Config *config, RunSummary *summary);

/* Get psi type label for display */
const char *analysis_psi_type_label(const Config *config);

/* Lookup known constant value by name
 * Returns true if found, false otherwise */
bool analysis_constant_value(const char *name, double *value);

#endif /* TRTS_ANALYSIS_UTILS_H */
