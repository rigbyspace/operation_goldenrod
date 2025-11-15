/* analysis_utils.c - Analysis Implementation
 *
 * In-memory observer-based analysis enforcing pure rational propagation.
 */

#include "analysis_utils.h"
#include "rational.h"
#include "simulate.h"
#include "state.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ARRAY_COUNT(arr) (sizeof(arr) / sizeof((arr)[0]))

/* Known mathematical constants for convergence detection */
typedef struct {
    const char *name;
    double value;
} KnownConstant;

static const KnownConstant KNOWN_CONSTANTS[] = {
    {"phi", 1.6180339887498948482},          /* Golden ratio */
    {"rho", 1.3247179572447458000},          /* Plastic constant */
    {"delta_s", 1.4655712318767680267},      /* Silver ratio */
    {"tribonacci", 1.8392867552141611326},   /* Tribonacci constant */
    {"plastic", 1.3247179572447458000},      /* Plastic (duplicate) */
    {"sqrt2", 1.4142135623730950488},        /* √2 */
    {"silver", 2.4142135623730950488}        /* Silver ratio (1+√2) */
};

/* In-memory observer context for collecting statistics */
typedef struct {
    RunSummary *summary;
    
    /* Welford algorithm state for ratio mean/variance */
    double ratio_mean;
    double ratio_m2;
    size_t ratio_count;
    double ratio_min;
    double ratio_max;
    
    size_t stack_sum;
    
    /* Delta and oscillation tracking */
    double previous_ratio;
    bool have_previous_ratio;
    double max_delta;
    size_t sign_changes;
    
    /* Best constant matching */
    double best_delta;
    size_t best_constant_index;
    
    /* Psi spacing (Welford) */
    size_t last_psi_index;
    bool have_last_psi;
    double spacing_mean;
    double spacing_m2;
    size_t spacing_count;
    
    /* Magnitude tracking */
    mpz_t max_mag_num;
    mpz_t max_mag_den;
    
    size_t last_tick_seen;
} InMemoryObserver;

/* Helper: convert mpq to double (for analysis only, never written back) */
static double mpq_to_double_snapshot(mpq_srcptr q) {
    return mpq_get_d(q);
}

/* Update max magnitude tracker */
static void update_max_mag(mpz_t max_mag, mpz_srcptr candidate) {
    mpz_t abs_candidate;
    mpz_init(abs_candidate);
    mpz_abs(abs_candidate, candidate);
    
    if (mpz_cmp(max_mag, abs_candidate) < 0) {
        mpz_set(max_mag, abs_candidate);
    }
    
    mpz_clear(abs_candidate);
}

/* Update stack summary string */
static void update_stack_summary(RunSummary *summary, size_t stack_sum) {
    if (summary->total_samples == 0U) {
        snprintf(summary->stack_summary, sizeof(summary->stack_summary), "avg=0.00 []");
        return;
    }
    
    summary->average_stack_depth = (double)stack_sum / (double)summary->total_samples;
    
    char buffer[128];
    int offset = snprintf(buffer, sizeof(buffer), "avg=%.2f [", summary->average_stack_depth);
    
    size_t histogram_limit = ARRAY_COUNT(summary->stack_histogram);
    for (size_t depth = 0; depth < histogram_limit; ++depth) {
        int written = snprintf(buffer + offset, sizeof(buffer) - (size_t)offset,
                              "%zu:%zu", depth, summary->stack_histogram[depth]);
        if (written < 0 || (size_t)written >= sizeof(buffer) - (size_t)offset) {
            break;
        }
        offset += written;
        
        if (depth + 1 < histogram_limit) {
            if ((size_t)offset + 1 < sizeof(buffer)) {
                buffer[offset++] = ',';
            }
        }
    }
    
    if ((size_t)offset + 1 < sizeof(buffer)) {
        buffer[offset++] = ']';
        buffer[offset] = '\0';
    } else {
        buffer[sizeof(buffer) - 1] = '\0';
    }
    
    strncpy(summary->stack_summary, buffer, sizeof(summary->stack_summary) - 1);
    summary->stack_summary[sizeof(summary->stack_summary) - 1] = '\0';
}

/* Determine pattern classification */
static void determine_pattern(RunSummary *summary, bool ratio_defined, bool divergent,
                              bool fixed_point, bool oscillating, size_t best_constant_index,
                              double best_delta) {
    if (!ratio_defined) {
        snprintf(summary->pattern, sizeof(summary->pattern), "null");
        snprintf(summary->classification, sizeof(summary->classification), "Null");
        return;
    }
    
    if (divergent) {
        snprintf(summary->pattern, sizeof(summary->pattern), "divergent");
        snprintf(summary->classification, sizeof(summary->classification), "Chaotic");
        return;
    }
    
    if (fixed_point) {
        snprintf(summary->pattern, sizeof(summary->pattern), "fixed point");
        snprintf(summary->classification, sizeof(summary->classification), "FixedPoint");
        return;
    }
    
    if (oscillating) {
        snprintf(summary->pattern, sizeof(summary->pattern), "oscillating");
        snprintf(summary->classification, sizeof(summary->classification), "Oscillating");
        return;
    }
    
    snprintf(summary->pattern, sizeof(summary->pattern), "stable");
    
    if (best_constant_index != SIZE_MAX && best_delta < 1e-4) {
        snprintf(summary->classification, sizeof(summary->classification), 
                "Convergent(%s)", KNOWN_CONSTANTS[best_constant_index].name);
    } else {
        snprintf(summary->classification, sizeof(summary->classification), "Stable");
    }
}

/* Observer callback implementation */
static void in_memory_observer(void *user_data, size_t tick, int microtick, char phase,
                               const TRTS_State *state, bool rho_event, bool psi_fired,
                               bool mu_zero, bool forced_emission) {
    InMemoryObserver *ctx = (InMemoryObserver *)user_data;
    RunSummary *summary = ctx->summary;
    
    /* Update tick counter */
    if (tick > ctx->last_tick_seen) {
        ctx->last_tick_seen = tick;
    }
    
    /* Update event counters */
    if (psi_fired) {
        summary->psi_events++;
        
        /* Compute psi spacing using linearized index */
        size_t current_index = (tick - 1U) * 11U + (size_t)microtick;
        if (ctx->have_last_psi) {
            double spacing = (double)(current_index - ctx->last_psi_index);
            ctx->spacing_count++;
            double delta = spacing - ctx->spacing_mean;
            ctx->spacing_mean += delta / (double)ctx->spacing_count;
            double delta2 = spacing - ctx->spacing_mean;
            ctx->spacing_m2 += delta * delta2;
        }
        ctx->last_psi_index = current_index;
        ctx->have_last_psi = true;
    }
    
    if (rho_event) {
        summary->rho_events++;
    }
    
    if (mu_zero) {
        summary->mu_zero_events++;
    }
    
    /* Update stack histogram */
    size_t stack_size = state->koppa_stack_size;
    if (stack_size >= ARRAY_COUNT(summary->stack_histogram)) {
        stack_size = ARRAY_COUNT(summary->stack_histogram) - 1U;
    }
    summary->stack_histogram[stack_size]++;
    ctx->stack_sum += stack_size;
    summary->total_samples++;
    
    /* Track max magnitude of numerators and denominators */
    update_max_mag(ctx->max_mag_num, mpq_numref(state->upsilon));
    update_max_mag(ctx->max_mag_den, mpq_denref(state->upsilon));
    update_max_mag(ctx->max_mag_num, mpq_numref(state->beta));
    update_max_mag(ctx->max_mag_den, mpq_denref(state->beta));
    
    /* Compute ratio if beta is non-zero */
    if (!rational_is_zero(&state->beta)) {
        Rational ratio_q;
        rational_init(&ratio_q);
        if (rational_div(&ratio_q, &state->upsilon, &state->beta)) {
            double snapshot = mpq_to_double_snapshot(ratio_q);
            
            summary->ratio_defined = true;
            summary->final_ratio_snapshot = snapshot;
            rational_set(summary->final_ratio, &ratio_q);
            
            /* Store textual representation */
            char *num_str = mpz_get_str(NULL, 10, mpq_numref(ratio_q));
            char *den_str = mpz_get_str(NULL, 10, mpq_denref(ratio_q));
            if (num_str && den_str) {
                snprintf(summary->final_ratio_str, sizeof(summary->final_ratio_str),
                        "%s/%s", num_str, den_str);
            }
            if (num_str) free(num_str);
            if (den_str) free(den_str);
            
            /* Welford update for ratio statistics */
            ctx->ratio_count++;
            if (ctx->ratio_count == 1U) {
                ctx->ratio_mean = snapshot;
                ctx->ratio_m2 = 0.0;
                ctx->ratio_min = snapshot;
                ctx->ratio_max = snapshot;
            } else {
                if (snapshot < ctx->ratio_min) ctx->ratio_min = snapshot;
                if (snapshot > ctx->ratio_max) ctx->ratio_max = snapshot;
                
                double delta = snapshot - ctx->ratio_mean;
                ctx->ratio_mean += delta / (double)ctx->ratio_count;
                double delta2 = snapshot - ctx->ratio_mean;
                ctx->ratio_m2 += delta * delta2;
            }
            
            /* Delta and oscillation tracking */
            if (ctx->have_previous_ratio) {
                double diff = fabs(snapshot - ctx->previous_ratio);
                if (diff > ctx->max_delta) {
                    ctx->max_delta = diff;
                }
                if ((snapshot > 0.0 && ctx->previous_ratio < 0.0) ||
                    (snapshot < 0.0 && ctx->previous_ratio > 0.0)) {
                    ctx->sign_changes++;
                }
            }
            ctx->previous_ratio = snapshot;
            ctx->have_previous_ratio = true;
            
            /* Check against known constants */
            for (size_t i = 0; i < ARRAY_COUNT(KNOWN_CONSTANTS); ++i) {
                double delta = fabs(snapshot - KNOWN_CONSTANTS[i].value);
                if (delta < ctx->best_delta) {
                    ctx->best_delta = delta;
                    ctx->best_constant_index = i;
                }
                if (delta < 1e-5 && summary->convergence_tick == 0U) {
                    summary->convergence_tick = tick;
                }
            }
        }
        rational_clear(&ratio_q);
    }
}

void run_summary_init(RunSummary *summary) {
    rational_init(summary->final_ratio);
    summary->ratio_defined = false;
    summary->final_ratio_str[0] = '\0';
    summary->closest_constant[0] = '\0';
    summary->closest_delta = 0.0;
    summary->convergence_tick = 0U;
    summary->pattern[0] = '\0';
    summary->classification[0] = '\0';
    summary->stack_summary[0] = '\0';
    summary->final_ratio_snapshot = 0.0;
    summary->total_samples = 0U;
    summary->total_ticks = 0U;
    summary->psi_events = 0U;
    summary->rho_events = 0U;
    summary->mu_zero_events = 0U;
    summary->psi_spacing_mean = 0.0;
    summary->psi_spacing_stddev = 0.0;
    summary->ratio_variance = 0.0;
    summary->ratio_range = 0.0;
    summary->ratio_mean = 0.0;
    summary->ratio_stddev = 0.0;
    summary->average_stack_depth = 0.0;
    for (size_t i = 0; i < ARRAY_COUNT(summary->stack_histogram); ++i) {
        summary->stack_histogram[i] = 0U;
    }
}

void run_summary_clear(RunSummary *summary) {
    rational_clear(summary->final_ratio);
}

void run_summary_copy(RunSummary *dest, const RunSummary *src) {
    if (dest == src) {
        return;
    }
    rational_set(dest->final_ratio, src->final_ratio);
    dest->ratio_defined = src->ratio_defined;
    memcpy(dest->final_ratio_str, src->final_ratio_str, sizeof(dest->final_ratio_str));
    memcpy(dest->closest_constant, src->closest_constant, sizeof(dest->closest_constant));
    dest->closest_delta = src->closest_delta;
    dest->convergence_tick = src->convergence_tick;
    memcpy(dest->pattern, src->pattern, sizeof(dest->pattern));
    memcpy(dest->classification, src->classification, sizeof(dest->classification));
    memcpy(dest->stack_summary, src->stack_summary, sizeof(dest->stack_summary));
    dest->final_ratio_snapshot = src->final_ratio_snapshot;
    dest->total_samples = src->total_samples;
    dest->total_ticks = src->total_ticks;
    dest->psi_events = src->psi_events;
    dest->rho_events = src->rho_events;
    dest->mu_zero_events = src->mu_zero_events;
    dest->psi_spacing_mean = src->psi_spacing_mean;
    dest->psi_spacing_stddev = src->psi_spacing_stddev;
    dest->ratio_variance = src->ratio_variance;
    dest->ratio_range = src->ratio_range;
    dest->ratio_mean = src->ratio_mean;
    dest->ratio_stddev = src->ratio_stddev;
    dest->average_stack_depth = src->average_stack_depth;
    memcpy(dest->stack_histogram, src->stack_histogram, sizeof(dest->stack_histogram));
}

bool analyze_latest_run(const Config *config, RunSummary *summary) {
    if (!config || !summary) {
        return false;
    }
    
    run_summary_init(summary);
    summary->total_ticks = 0U;
    
    /* Initialize observer context */
    InMemoryObserver ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.summary = summary;
    ctx.best_delta = INFINITY;
    ctx.best_constant_index = SIZE_MAX;
    mpz_init_set_ui(ctx.max_mag_num, 0UL);
    mpz_init_set_ui(ctx.max_mag_den, 0UL);
    
    /* Run simulation with observer */
    simulate_stream(config, in_memory_observer, &ctx);
    
    /* Finalize statistics */
    summary->total_ticks = ctx.last_tick_seen;
    
    if (ctx.ratio_count > 0U) {
        summary->ratio_mean = ctx.ratio_mean;
        if (ctx.ratio_count > 1U) {
            summary->ratio_variance = ctx.ratio_m2 / (double)(ctx.ratio_count - 1U);
            summary->ratio_stddev = sqrt(summary->ratio_variance);
        }
        summary->ratio_range = ctx.ratio_max - ctx.ratio_min;
    }
    
    if (ctx.spacing_count > 1U) {
        summary->psi_spacing_mean = ctx.spacing_mean;
        summary->psi_spacing_stddev = sqrt(ctx.spacing_m2 / (double)(ctx.spacing_count - 1U));
    } else if (ctx.spacing_count == 1U) {
        summary->psi_spacing_mean = ctx.spacing_mean;
    }
    
    update_stack_summary(summary, ctx.stack_sum);
    
    /* Store closest constant */
    if (ctx.best_constant_index != SIZE_MAX) {
        strncpy(summary->closest_constant, KNOWN_CONSTANTS[ctx.best_constant_index].name,
               sizeof(summary->closest_constant) - 1);
        summary->closest_constant[sizeof(summary->closest_constant) - 1] = '\0';
        summary->closest_delta = ctx.best_delta;
    } else {
        strncpy(summary->closest_constant, "None", sizeof(summary->closest_constant) - 1);
        summary->closest_constant[sizeof(summary->closest_constant) - 1] = '\0';
        summary->closest_delta = INFINITY;
    }
    
    /* Determine pattern classification */
    mpz_t divergence_threshold;
    mpz_init_set_ui(divergence_threshold, 1000000000UL);
    
    bool divergent = summary->ratio_defined &&
                    (summary->ratio_range > 1.0e6 ||
                     mpz_cmp(ctx.max_mag_num, divergence_threshold) > 0 ||
                     mpz_cmp(ctx.max_mag_den, divergence_threshold) > 0);
    
    bool fixed_point = summary->ratio_defined && 
                      summary->ratio_range < 1.0e-9 && 
                      ctx.max_delta < 1.0e-12;
    
    bool oscillating = summary->ratio_defined && !divergent && !fixed_point &&
                      summary->ratio_range < 100.0 && 
                      ctx.sign_changes > ctx.ratio_count / 3U;
    
    determine_pattern(summary, summary->ratio_defined, divergent, fixed_point, 
                     oscillating, ctx.best_constant_index, ctx.best_delta);
    
    mpz_clear(divergence_threshold);
    mpz_clear(ctx.max_mag_num);
    mpz_clear(ctx.max_mag_den);
    
    return true;
}

bool simulate_and_analyze(const Config *config, RunSummary *summary) {
    return analyze_latest_run(config, summary);
}

const char *analysis_psi_type_label(const Config *config) {
    return config->triple_psi_mode ? "3-way" : "2-way";
}

bool analysis_constant_value(const char *name, double *value) {
    if (!name || !value) {
        return false;
    }
    for (size_t i = 0; i < ARRAY_COUNT(KNOWN_CONSTANTS); ++i) {
        if (strcmp(name, KNOWN_CONSTANTS[i].name) == 0) {
            *value = KNOWN_CONSTANTS[i].value;
            return true;
        }
    }
    return false;
}
