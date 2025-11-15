/* analysis_utils.c - Analysis Implementation
 *
 * In-memory observer-based analysis enforcing pure rational propagation.
 */
#include "stdint.h"
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
    {"silver", 2.4142135623730950488},       /* Silver ratio alternate */
    {"e", 2.7182818284590452354},            /* Euler's number */
    {"pi", 3.1415926535897932384}             /* Pi */
};

/* Helper: Update maximum magnitude tracking */
static void update_max_mag(mpz_t current_max, const mpz_t new_value) {
    mpz_t abs_val;
    mpz_init(abs_val);
    mpz_abs(abs_val, new_value);
    
    if (mpz_cmp(abs_val, current_max) > 0) {
        mpz_set(current_max, abs_val);
    }
    
    mpz_clear(abs_val);
}

/* Helper: Performs one step of Welford's algorithm for online variance */
static void welford_update(double *count, double *mean, double *m2, double new_value) {
    double delta, delta2;
    
    (*count)++;
    delta = new_value - *mean;
    *mean += delta / *count;
    delta2 = new_value - *mean;
    *m2 += delta * delta2;
}

/* Helper: Calculates variance and standard deviation from Welford state */
static void welford_finalize(double count, double m2, double *variance, double *stddev) {
    if (count > 1) {
        *variance = m2 / (count - 1.0);
        *stddev = sqrt(*variance);
    } else {
        *variance = 0.0;
        *stddev = 0.0;
    }
}

/* Helper: Determines the pattern classification */
static void determine_pattern(RunSummary *summary, bool defined, bool divergent, 
                              bool fixed_point, bool oscillating, int best_idx, double best_delta) {
    if (!defined) {
        strncpy(summary->pattern, "undefined", sizeof(summary->pattern));
        snprintf(summary->classification, sizeof(summary->classification), 
                 "Ratio undefined (division by zero encountered)");
        return;
    }
    
    if (divergent) {
        strncpy(summary->pattern, "divergent", sizeof(summary->pattern));
        snprintf(summary->classification, sizeof(summary->classification), 
                 "Magnitude growing rapidly (max component > 1e6 or 1e10)");
        return;
    }
    
    if (fixed_point) {
        strncpy(summary->pattern, "fixed point", sizeof(summary->pattern));
        snprintf(summary->classification, sizeof(summary->classification), 
                 "Converged to a stable value (ratio range < 1e-9)");
    } else if (oscillating) {
        strncpy(summary->pattern, "oscillating", sizeof(summary->pattern));
        snprintf(summary->classification, sizeof(summary->classification), 
                 "Oscillating or cyclic behavior (high sign flip frequency)");
    } else {
        strncpy(summary->pattern, "stable", sizeof(summary->pattern));
        snprintf(summary->classification, sizeof(summary->classification), 
                 "Stable and non-divergent");
    }
    
    if (best_idx >= 0 && best_delta < 1.0e-5) {
        const char *constant_name = KNOWN_CONSTANTS[best_idx].name;
        snprintf(summary->classification, sizeof(summary->classification), 
                 "%s, potentially converging to %s", summary->classification, constant_name);
        strncpy(summary->closest_constant, constant_name, sizeof(summary->closest_constant));
        summary->closest_delta = best_delta;
    }
}

/* Context for in-memory observation */
typedef struct {
    RunSummary *summary;
    size_t tick_count;
    mpz_t max_mag_num;
    mpz_t max_mag_den;
    
    /* Welford state for ratio variance */
    double ratio_count;
    double ratio_mean;
    double ratio_m2;
    double ratio_min;
    double ratio_max;
    
    /* Welford state for psi spacing */
    double psi_spacing_count;
    double psi_spacing_mean;
    double psi_spacing_m2;
    int last_psi_mt;
    
    /* Convergence tracking */
    int best_constant_index;
    double best_delta;
    size_t convergence_tick;
    
    /* Other metrics */
    size_t psi_fire_count;
    size_t psi_triple_count;
    int sign_changes;
    double max_delta;
} AnalysisContext;

/* Observer callback used during simulation */
bool in_memory_observer(void *user_data, size_t tick, int microtick, char phase,
                        const TRTS_State *state, const Config *config,
                        bool psi_fired, bool koppa_sampled,
                        bool mu_zero, bool forced_emission) {

    AnalysisContext *ctx = (AnalysisContext *)user_data;
    RunSummary *summary = ctx->summary;
    
    /* --- Phase 1: Engine Step/Magnitude Tracking (Microticks 1-10) --- */
    if (microtick >= 1 && microtick <= 10) {
        ctx->tick_count++;
        summary->engine_step_count++;
        
        /* Magnitude tracking for divergence detection */
        // FIX: Access struct members directly instead of GMP macros
        update_max_mag(ctx->max_mag_num, state->upsilon.num); 
        update_max_mag(ctx->max_mag_den, state->upsilon.den); 
        update_max_mag(ctx->max_mag_num, state->beta.num);    
        update_max_mag(ctx->max_mag_den, state->beta.den);    
        
        /* Track deltas */
        double delta_u = mpz_get_d(state->delta_upsilon.num) / mpz_get_d(state->delta_upsilon.den);
        double delta_b = mpz_get_d(state->delta_beta.num) / mpz_get_d(state->delta_beta.den);
        double delta_mag = fabs(delta_u) + fabs(delta_b);
        if (delta_mag > ctx->max_delta) {
            ctx->max_delta = delta_mag;
        }
    }
    
    /* --- Phase 2: Ratio Sampling (Microtick 10 only) --- */
    if (microtick == 10 && !rational_denominator_zero(&state->beta)) {
        
        Rational ratio_q;
        rational_init(&ratio_q);
        
        if (rational_div(&ratio_q, &state->upsilon, &state->beta)) {
            summary->ratio_defined = true;
            
            // FIX: Direct division using mpz_get_d instead of GMP-specific function
            double snapshot = mpz_get_d(ratio_q.num) / mpz_get_d(ratio_q.den); 
            
            /* Welford update */
            welford_update(&ctx->ratio_count, &ctx->ratio_mean, &ctx->ratio_m2, snapshot);
            if (snapshot < ctx->ratio_min) ctx->ratio_min = snapshot;
            if (snapshot > ctx->ratio_max) ctx->ratio_max = snapshot;
            
            /* Check for nearest constant */
            for (int i = 0; i < ARRAY_COUNT(KNOWN_CONSTANTS); ++i) {
                double delta = fabs(snapshot - KNOWN_CONSTANTS[i].value);
                if (delta < ctx->best_delta) {
                    ctx->best_delta = delta;
                    ctx->best_constant_index = i;
                    if (delta < 1.0e-5 && ctx->convergence_tick == 0) {
                        ctx->convergence_tick = tick;
                    }
                }
            }
            
            /* Store last ratio as string */
            // FIX: Pass address (&) of the Rational struct member
            rational_set(&summary->final_ratio, &ratio_q); 
            
            // FIX: Access struct members directly instead of GMP macros
            char *num_str = mpz_get_str(NULL, 10, ratio_q.num); 
            char *den_str = mpz_get_str(NULL, 10, ratio_q.den); 
            snprintf(summary->final_ratio_str, sizeof(summary->final_ratio_str),
                     "%s/%s", num_str, den_str);
            free(num_str);
            free(den_str);

            /* Track sign changes */
            if (ctx->ratio_count > 1.0) {
                if ((snapshot > 0.0 && summary->ratio_mean < 0.0) || 
                    (snapshot < 0.0 && summary->ratio_mean > 0.0)) {
                    ctx->sign_changes++;
                }
            }
        } else {
            summary->ratio_defined = false;
        }
        
        rational_clear(&ratio_q);
    }
    
    /* --- Phase 3: Koopa Sampling (Microtick 11 only) --- */
    if (microtick == 11 && koppa_sampled) {
        summary->koppa_sample_count++;
        /* ... (Koopa stack summary logic goes here) ... */
        
        /* For simplicity, just store the size */
        snprintf(summary->stack_summary, sizeof(summary->stack_summary), 
                 "Depth: %zu", state->koppa_stack_size);
        if (state->koppa_stack_size > summary->stack_max_depth) {
            summary->stack_max_depth = state->koppa_stack_size;
        }
    }
    
    /* --- Phase 4: Psi Firing Tracking (Any Microtick) --- */
    if (psi_fired) {
        summary->psi_fire_count++;
        if (state->psi_triple_recent) {
            summary->psi_triple_count++;
        }
        
        /* Welford update for spacing */
        int spacing = microtick - ctx->last_psi_mt;
        if (ctx->last_psi_mt != 0) {
            welford_update(&ctx->psi_spacing_count, &ctx->psi_spacing_mean, &ctx->psi_spacing_m2, (double)spacing);
        }
        ctx->last_psi_mt = microtick;
    }
    
    return true; /* Continue simulation */
}

void run_summary_init(RunSummary *summary) {
    if (!summary) return;
    
    // FIX: Pass address (&) of the Rational struct member
    rational_init(&summary->final_ratio); 
    
    summary->ratio_defined = false;
    summary->final_ratio_str[0] = '\0';
    
    summary->closest_constant[0] = '\0';
    summary->closest_delta = 0.0;
    summary->convergence_tick = 0;
    
    summary->pattern[0] = '\0';
    summary->classification[0] = '\0';
    
    summary->stack_summary[0] = '\0';
    summary->stack_max_depth = 0;
    
    summary->engine_step_count = 0;
    summary->koppa_sample_count = 0;
    
    summary->psi_fire_count = 0;
    summary->psi_triple_count = 0;
    summary->psi_spacing_mean = 0.0;
    summary->psi_spacing_stddev = 0.0;
    
    summary->ratio_variance = 0.0;
    summary->ratio_range = 0.0;
    summary->ratio_mean = 0.0;
    summary->ratio_stddev = 0.0;
}

void run_summary_clear(RunSummary *summary) {
    if (!summary) return;
    
    // FIX: Pass address (&) of the Rational struct member
    rational_clear(&summary->final_ratio); 
}

void run_summary_copy(RunSummary *dest, const RunSummary *src) {
    if (!dest || !src) return;
    
    // FIX: Pass addresses (&) of both Rational struct members
    rational_set(&dest->final_ratio, &src->final_ratio); 
    
    dest->ratio_defined = src->ratio_defined;
    strncpy(dest->final_ratio_str, src->final_ratio_str, sizeof(dest->final_ratio_str));
    
    dest->closest_delta = src->closest_delta;
    dest->convergence_tick = src->convergence_tick;
    strncpy(dest->closest_constant, src->closest_constant, sizeof(dest->closest_constant));
    
    strncpy(dest->pattern, src->pattern, sizeof(dest->pattern));
    strncpy(dest->classification, src->classification, sizeof(dest->classification));
    
    strncpy(dest->stack_summary, src->stack_summary, sizeof(dest->stack_summary));
    dest->stack_max_depth = src->stack_max_depth;
    
    dest->engine_step_count = src->engine_step_count;
    dest->koppa_sample_count = src->koppa_sample_count;
    
    dest->psi_fire_count = src->psi_fire_count;
    dest->psi_triple_count = src->psi_triple_count;
    dest->psi_spacing_mean = src->psi_spacing_mean;
    dest->psi_spacing_stddev = src->psi_spacing_stddev;
    
    dest->ratio_variance = src->ratio_variance;
    dest->ratio_range = src->ratio_range;
    dest->ratio_mean = src->ratio_mean;
    dest->ratio_stddev = src->ratio_stddev;
}

bool analyze_latest_run(const Config *config, RunSummary *summary) {
    if (!config || !summary) {
        return false;
    }
    
    AnalysisContext ctx = {0};
    ctx.summary = summary;
    ctx.ratio_min = INFINITY;
    ctx.ratio_max = -INFINITY;
    ctx.best_delta = INFINITY;
    ctx.best_constant_index = -1;
    mpz_init(ctx.max_mag_num);
    mpz_init(ctx.max_mag_den);
    
    /* Reset summary fields that are cumulative */
    summary->engine_step_count = 0;
    summary->koppa_sample_count = 0;
    summary->psi_fire_count = 0;
    summary->psi_triple_count = 0;
    summary->stack_max_depth = 0;
    
    /* Re-initialize GMP rational in case it was cleared */
    rational_init(&summary->final_ratio);

    /* Run simulation with observer */
    simulate_stream(config, in_memory_observer, &ctx);

    /* Finalize statistics */
    welford_finalize(ctx.psi_spacing_count, ctx.psi_spacing_m2, 
                     &summary->psi_spacing_mean, &summary->psi_spacing_stddev);

    welford_finalize(ctx.ratio_count, ctx.ratio_m2, 
                     &summary->ratio_variance, &summary->ratio_stddev);
    
    summary->ratio_mean = ctx.ratio_mean;
    if (ctx.ratio_count > 0) {
        summary->ratio_range = ctx.ratio_max - ctx.ratio_min;
    } else {
        summary->ratio_range = 0.0;
    }
    
    /* Check for divergence */
    mpz_t divergence_threshold;
    mpz_init_set_ui(divergence_threshold, 10000000000UL); // 10^10
    
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
    if (!name || !value) return false;
    
    for (int i = 0; i < ARRAY_COUNT(KNOWN_CONSTANTS); ++i) {
        if (strcmp(name, KNOWN_CONSTANTS[i].name) == 0) {
            *value = KNOWN_CONSTANTS[i].value;
            return true;
        }
    }
    return false;
}
