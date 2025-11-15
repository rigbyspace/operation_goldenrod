// analysis_utils.h
#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>

#include "config.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    mpq_t final_ratio;
    bool ratio_defined;
    char final_ratio_str[128];
    char closest_constant[32];
    double closest_delta;
    size_t convergence_tick;
    char pattern[32];
    char classification[64];
    char stack_summary[128];
    double final_ratio_snapshot;
    size_t total_samples;
    size_t total_ticks;
    size_t psi_events;
    size_t rho_events;
    size_t mu_zero_events;
    double psi_spacing_mean;
    double psi_spacing_stddev;
    double ratio_variance;
    double ratio_range;
    double ratio_mean;
    double ratio_stddev;
    size_t stack_histogram[8];
    double average_stack_depth;
} RunSummary;

void run_summary_init(RunSummary *summary);
void run_summary_clear(RunSummary *summary);
void run_summary_copy(RunSummary *dest, const RunSummary *src);

bool analyze_latest_run(const Config *config, RunSummary *summary);
bool simulate_and_analyze(const Config *config, RunSummary *summary);

const char *analysis_psi_type_label(const Config *config);
bool analysis_constant_value(const char *name, double *value);

bool simulate_and_analyze(const Config *config, RunSummary *summary);


#ifdef __cplusplus
}
#endif

#endif // ANALYSIS_UTILS_H

