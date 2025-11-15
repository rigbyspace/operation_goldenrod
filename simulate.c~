// simulate.c
/*
 * Complete simulation loop with all helper functions implemented.
 */

#include "simulate.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>

#include "engine.h"
#include "koppa.h"
#include "psi.h"
#include "rational.h"

/* ===========================================================
   PRIME AND PATTERN CHECK LOGIC
   =========================================================== */

static bool mpz_is_prime_signed(mpz_srcptr value) {
    mpz_t magnitude;
    mpz_init(magnitude);
    mpz_abs(magnitude, value);
    
    bool is_prime = false;
    if (mpz_cmp_ui(magnitude, 2UL) >= 0) {
        is_prime = mpz_probab_prime_p(magnitude, 25) > 0;
    }
    mpz_clear(magnitude);
    return is_prime;
}

static bool mpz_is_square(mpz_srcptr value) {
    if (mpz_sgn(value) < 0) {
        return false;
    }
    mpz_t root;
    mpz_init(root);
    mpz_sqrt(root, value);
    mpz_mul(root, root, root);
    bool result = (mpz_cmp(root, value) == 0);
    mpz_clear(root);
    return result;
}

static bool mpz_is_fibonacci(mpz_srcptr value) {
    if (mpz_cmp_ui(value, 0UL) < 0) {
        return false;
    }
    if (mpz_cmp_ui(value, 1UL) <= 0) {
        return true;
    }
    
    mpz_t test1, test2, temp;
    mpz_init(test1);
    mpz_init(test2);
    mpz_init(temp);
    
    // Test: n is Fibonacci iff one of 5*n^2 + 4 or 5*n^2 - 4 is a perfect square
    mpz_mul(temp, value, value);
    mpz_mul_ui(temp, temp, 5UL);
    
    mpz_add_ui(test1, temp, 4UL);
    mpz_sub_ui(test2, temp, 4UL);
    
    bool result = mpz_is_square(test1) || mpz_is_square(test2);
    
    mpz_clear(test1);
    mpz_clear(test2);
    mpz_clear(temp);
    return result;
}

static bool mpz_is_perfect_power(mpz_srcptr value) {
    if (mpz_cmp_ui(value, 1UL) <= 0) {
        return false;
    }
    return mpz_perfect_power_p(value) != 0;
}

static bool mpq_has_pattern_component(const Config *config, mpq_srcptr value,
                                       bool check_num, bool check_den) {
    bool found = false;
    
    if (check_num) {
        mpz_srcptr num = mpq_numref(value);
        if (mpz_is_prime_signed(num)) {
            found = true;
        }
        if (config->enable_twin_prime_trigger && mpz_is_prime_signed(num)) {
            mpz_t temp;
            mpz_init(temp);
            mpz_add_ui(temp, num, 2UL);
            if (mpz_is_prime_signed(temp)) {
                found = true;
            }
            mpz_sub_ui(temp, num, 2UL);
            if (mpz_is_prime_signed(temp)) {
                found = true;
            }
            mpz_clear(temp);
        }
        if (config->enable_fibonacci_trigger) {
            mpz_t abs_num;
            mpz_init(abs_num);
            mpz_abs(abs_num, num);
            if (mpz_is_fibonacci(abs_num)) {
                found = true;
            }
            mpz_clear(abs_num);
        }
        if (config->enable_perfect_power_trigger) {
            mpz_t abs_num;
            mpz_init(abs_num);
            mpz_abs(abs_num, num);
            if (mpz_is_perfect_power(abs_num)) {
                found = true;
            }
            mpz_clear(abs_num);
        }
    }
    
    if (check_den) {
        mpz_srcptr den = mpq_denref(value);
        if (mpz_is_prime_signed(den)) {
            found = true;
        }
        if (config->enable_fibonacci_trigger && mpz_is_fibonacci(den)) {
            found = true;
        }
        if (config->enable_perfect_power_trigger && mpz_is_perfect_power(den)) {
            found = true;
        }
    }
    
    return found;
}

/* ===========================================================
   RATIO TRIGGERS AND PSI CONDITIONS
   =========================================================== */

static void ratio_bounds(RatioTriggerMode mode, mpq_t lower, mpq_t upper) {
    switch (mode) {
    case RATIO_TRIGGER_GOLDEN:
        rational_set_si(lower, 3, 2);
        rational_set_si(upper, 17, 10);
        break;
    case RATIO_TRIGGER_SQRT2:
        rational_set_si(lower, 13, 10);
        rational_set_si(upper, 3, 2);
        break;
    case RATIO_TRIGGER_PLASTIC:
        rational_set_si(lower, 6, 5);
        rational_set_si(upper, 7, 5);
        break;
    case RATIO_TRIGGER_NONE:
    case RATIO_TRIGGER_CUSTOM:
    default:
        rational_set_si(lower, 0, 1);
        rational_set_si(upper, 0, 1);
        break;
    }
}

static bool ratio_in_range(const Config *config, const TRTS_State *state) {
    if (config->ratio_trigger_mode == RATIO_TRIGGER_NONE) {
        return false;
    }
    if (rational_is_zero(state->beta)) {
        return false;
    }
    mpq_t ratio;
    rational_init(ratio);
    rational_div(ratio, state->upsilon, state->beta);
    bool in_range = false;
    if (config->ratio_trigger_mode == RATIO_TRIGGER_CUSTOM && config->enable_ratio_custom_range) {
        if (mpq_cmp(ratio, config->ratio_custom_lower) > 0 && 
            mpq_cmp(ratio, config->ratio_custom_upper) < 0) {
            in_range = true;
        }
    } else {
        mpq_t lower, upper;
        rational_init(lower);
        rational_init(upper);
        ratio_bounds(config->ratio_trigger_mode, lower, upper);
        if (mpq_cmp(ratio, lower) > 0 && mpq_cmp(ratio, upper) < 0) {
            in_range = true;
        }
        rational_clear(lower);
        rational_clear(upper);
    }
    rational_clear(ratio);
    return in_range;
}

static bool ratio_threshold_outside(const Config *config, const TRTS_State *state) {
    if (!config->enable_ratio_threshold_psi) {
        return false;
    }
    if (rational_is_zero(state->beta)) {
        return false;
    }
    mpq_t ratio;
    rational_init(ratio);
    rational_div(ratio, state->upsilon, state->beta);
    double ratio_snapshot = mpq_get_d(ratio);
    rational_clear(ratio);
    double magnitude = ratio_snapshot >= 0.0 ? ratio_snapshot : -ratio_snapshot;
    return magnitude < 0.5 || magnitude > 2.0;
}

static bool should_fire_psi(const Config *config, const TRTS_State *state,
                            bool is_memory_step, bool allow_stack) {
    if (!is_memory_step || !allow_stack) {
        return false;
    }
    switch (config->psi_mode) {
    case PSI_MODE_MSTEP:
        return true;
    case PSI_MODE_RHO_ONLY:
        return state->rho_pending;
    case PSI_MODE_MSTEP_RHO:
        return true;
    case PSI_MODE_INHIBIT_RHO:
        return !state->rho_pending;
    }
    return false;
}

static bool stack_allows_psi(const Config *config, const TRTS_State *state) {
    if (!config->enable_stack_depth_modes) {
        return true;
    }
    return state->koppa_stack_size == 2 || state->koppa_stack_size == 4;
}

/* ===========================================================
   SIMULATION OUTPUT HANDLING
   =========================================================== */

typedef struct {
    FILE *events_file;
    FILE *values_file;
} SimulationOutputs;

static void log_event(FILE *events_file, size_t tick, int microtick, char phase,
                       bool rho_event, bool psi_fired, bool mu_zero, bool forced_emission,
                       const TRTS_State *state) {
    fprintf(events_file,
            "%zu,%d,%c,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",
            tick, microtick, phase,
            rho_event ? 1 : 0, psi_fired ? 1 : 0, mu_zero ? 1 : 0,
            forced_emission ? 1 : 0,
            state->ratio_triggered_recent ? 1 : 0, state->psi_triple_recent ? 1 : 0,
            state->dual_engine_last_step ? 1 : 0, state->koppa_sample_index,
            state->ratio_threshold_recent ? 1 : 0, state->psi_strength_applied ? 1 : 0,
            state->sign_flip_polarity ? 1 : 0);
}

static void log_values(FILE *values_file, size_t tick, int microtick,
                        const TRTS_State *state) {
    gmp_fprintf(values_file,
        "%zu,%d,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%zu,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",
        tick, microtick,
        mpq_numref(state->upsilon), mpq_denref(state->upsilon),
        mpq_numref(state->beta), mpq_denref(state->beta),
        mpq_numref(state->koppa), mpq_denref(state->koppa),
        mpq_numref(state->koppa_sample), mpq_denref(state->koppa_sample),
        mpq_numref(state->previous_upsilon), mpq_denref(state->previous_upsilon),
        mpq_numref(state->previous_beta), mpq_denref(state->previous_beta),
        mpq_numref(state->koppa_stack[0]), mpq_denref(state->koppa_stack[0]),
        mpq_numref(state->koppa_stack[1]), mpq_denref(state->koppa_stack[1]),
        mpq_numref(state->koppa_stack[2]), mpq_denref(state->koppa_stack[2]),
        mpq_numref(state->koppa_stack[3]), mpq_denref(state->koppa_stack[3]),
        state->koppa_stack_size,
        mpq_numref(state->delta_upsilon), mpq_denref(state->delta_upsilon),
        mpq_numref(state->delta_beta), mpq_denref(state->delta_beta),
        mpq_numref(state->triangle_phi_over_epsilon), mpq_denref(state->triangle_phi_over_epsilon),
        mpq_numref(state->triangle_prev_over_phi), mpq_denref(state->triangle_prev_over_phi),
        mpq_numref(state->triangle_epsilon_over_prev), mpq_denref(state->triangle_epsilon_over_prev));
}

static void emit_outputs(const SimulationOutputs *outputs, size_t tick, int microtick,
                          char phase, bool rho_event, bool psi_fired, bool mu_zero,
                          bool forced_emission, const TRTS_State *state,
                          SimulateObserver observer, void *user_data) {
    if (outputs) {
        if (outputs->events_file) {
            log_event(outputs->events_file, tick, microtick, phase,
                      rho_event, psi_fired, mu_zero, forced_emission, state);
        }
        if (outputs->values_file) {
            log_values(outputs->values_file, tick, microtick, state);
        }
    }
    if (observer) {
        observer(user_data, tick, microtick, phase, state,
                 rho_event, psi_fired, mu_zero, forced_emission);
    }
}

/* ===========================================================
   CORE SIMULATION LOOP
   =========================================================== */

static void run_simulation(const Config *config, const SimulationOutputs *outputs,
                            SimulateObserver observer, void *user_data) {
    TRTS_State state;
    state_init(&state);
    state_reset(&state, config);
    
    for (size_t tick = 1; tick <= config->ticks; ++tick) {
        state.tick = tick;
        for (int microtick = 1; microtick <= 11; ++microtick) {
            char phase;
            switch (microtick) {
            case 1:
            case 4:
            case 7:
            case 10:
                phase = 'E';
                break;
            case 2:
            case 5:
            case 8:
            case 11:
                phase = 'M';
                break;
            default:
                phase = 'R';
                break;
            }
            
            bool rho_event = false;
            bool psi_fired = false;
            bool mu_zero = false;
            bool forced_emission = false;
            
            // Clear per-microtick flags
            state.ratio_triggered_recent = false;
            state.psi_triple_recent = false;
            state.dual_engine_last_step = false;
            state.koppa_sample_index = -1;
            rational_set(state.koppa_sample, state.koppa);
            state.ratio_threshold_recent = false;
            state.psi_strength_applied = false;
            
            switch (phase) {
            case 'E': {
                // Epsilon phase: compute epsilon and run engine step
                rational_set(state.epsilon, state.upsilon);
                bool engine_ok = engine_step(config, &state, microtick);
                (void)engine_ok;

                // Check for pattern components if enabled
                if (config->prime_target == PRIME_ON_NEW_UPSILON) {
                    if (mpq_has_pattern_component(config, state.upsilon, true, false)) {
                        state.rho_pending = true;
                        rho_event = true;
                    }
                }

                // Microtick 10 behavior
                forced_emission = (microtick == 10);
                if (microtick == 10 && config->mt10_behavior == MT10_FORCED_PSI) {
                    state.rho_pending = true;
                    rho_event = true;
                }
                break;
            }
            case 'M': {
                // Memory phase: decide whether to fire ψ
                mu_zero = rational_is_zero(state.beta);
                
                // Check memory value for patterns
                if (config->prime_target == PRIME_ON_MEMORY) {
                    if (mpq_has_pattern_component(config, state.beta, true, true)) {
                        state.rho_pending = true;
                        rho_event = true;
                    }
                }
                
                bool allow_stack = stack_allows_psi(config, &state);
                bool request_psi = should_fire_psi(config, &state, true, allow_stack);
                bool ratio_triggered = ratio_in_range(config, &state);
                if (ratio_triggered) {
                    request_psi = true;
                }
                bool ratio_threshold = ratio_threshold_outside(config, &state);
                if (ratio_threshold) {
                    request_psi = true;
                    state.ratio_threshold_recent = true;
                }
                if (request_psi && allow_stack) {
                    psi_fired = psi_transform(config, &state);
                } else {
                    state.psi_recent = false;
                }
                state.ratio_triggered_recent = ratio_triggered;
                
                // Accrue κ and reset ρ latch for the next microtick
                koppa_accrue(config, &state, psi_fired, true, microtick);
                state.rho_latched = false;
                break;
            }
            case 'R': {
                // Reset phase: accrue κ without psi
                koppa_accrue(config, &state, false, false, microtick);
                state.psi_recent = false;
                state.rho_latched = false;
                break;
            }
            }
            
            // Emit outputs (log or observer callback)
            emit_outputs(outputs, tick, microtick, phase, rho_event,
                         psi_fired, mu_zero, forced_emission, &state,
                         observer, user_data);
        }
    }
    state_clear(&state);
}

/* ===========================================================
   ENTRY POINTS
   =========================================================== */

void simulate(const Config *config) {
    FILE *events_file = fopen("events.csv", "w");
    if (!events_file) {
        perror("events.csv");
        return;
    }
    FILE *values_file = fopen("values.csv", "w");
    if (!values_file) {
        perror("values.csv");
        fclose(events_file);
        return;
    }
    
    // Write CSV headers
    fprintf(events_file,
            "tick,mt,phase,rho_event,psi_fired,mu_zero,forced_emission,"
            "ratio_triggered,triple_psi,dual_engine,koppa_sample_index,"
            "ratio_threshold,psi_strength,sign_flip\n");
    fprintf(values_file,
            "tick,mt,upsilon_num,upsilon_den,beta_num,beta_den,koppa_num,koppa_den,"
            "koppa_sample_num,koppa_sample_den,prev_upsilon_num,prev_upsilon_den,"
            "prev_beta_num,prev_beta_den,koppa_stack0_num,koppa_stack0_den,"
            "koppa_stack1_num,koppa_stack1_den,koppa_stack2_num,koppa_stack2_den,"
            "koppa_stack3_num,koppa_stack3_den,koppa_stack_size,delta_upsilon_num,"
            "delta_upsilon_den,delta_beta_num,delta_beta_den,triangle_phi_over_epsilon_num,"
            "triangle_phi_over_epsilon_den,triangle_prev_over_phi_num,"
            "triangle_prev_over_phi_den,triangle_epsilon_over_prev_num,"
            "triangle_epsilon_over_prev_den\n");
    
    SimulationOutputs outputs = {events_file, values_file};
    run_simulation(config, &outputs, NULL, NULL);
    
    fclose(events_file);
    fclose(values_file);
}

void simulate_stream(const Config *config, SimulateObserver observer, void *user_data) {
    run_simulation(config, NULL, observer, user_data);
}

