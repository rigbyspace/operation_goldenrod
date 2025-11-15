// self_refine.c
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "simulate.h"
#include "analysis_utils.h"
#include "config.h"

#define ARRAY_COUNT(arr) (sizeof(arr) / sizeof((arr)[0]))

typedef struct {
    Config config;
    RunSummary summary;
    double score;
    bool evaluated;
} Candidate;

typedef struct {
    unsigned generations;
    unsigned population;
    unsigned elite;
    unsigned seed;
    char summary_dir[256];
    bool use_summary_dir;
} EvolutionOptions;

static void candidate_init(Candidate *candidate) {
    config_init(&candidate->config);
    run_summary_init(&candidate->summary);
    candidate->score = 0.0;
    candidate->evaluated = false;
}

static void candidate_clear(Candidate *candidate) {
    config_clear(&candidate->config);
    run_summary_clear(&candidate->summary);
}

static void config_clone(Config *dest, const Config *src) {
    config_clear(dest);
    config_init(dest);

    dest->ticks = src->ticks;
    dest->engine_mode = src->engine_mode;
    dest->engine_upsilon = src->engine_upsilon;
    dest->engine_beta = src->engine_beta;
    dest->psi_mode = src->psi_mode;
    dest->koppa_mode = src->koppa_mode;
    dest->koppa_trigger = src->koppa_trigger;
    dest->prime_target = src->prime_target;
    dest->mt10_behavior = src->mt10_behavior;

    rational_set(dest->initial_koppa, src->initial_koppa);
    rational_set(dest->psi_alpha, src->psi_alpha);
    rational_set(dest->psi_beta, src->psi_beta);
    rational_set(dest->psi_gamma, src->psi_gamma);
    rational_set(dest->psi_delta, src->psi_delta);
    rational_set(dest->psi_epsilon, src->psi_epsilon);
    rational_set(dest->koppa_alpha, src->koppa_alpha);
    rational_set(dest->koppa_beta, src->koppa_beta);
    rational_set(dest->koppa_gamma, src->koppa_gamma);
    rational_set(dest->koppa_delta, src->koppa_delta);
    rational_set(dest->koppa_epsilon, src->koppa_epsilon);
    rational_set(dest->koppa_zeta, src->koppa_zeta);
    rational_set(dest->koppa_eta, src->koppa_eta);
    rational_set(dest->koppa_theta, src->koppa_theta);
    rational_set(dest->koppa_iota, src->koppa_iota);
    rational_set(dest->koppa_kappa, src->koppa_kappa);
}

static double random_double(double min_value, double max_value) {
    double t = (double)rand() / (double)RAND_MAX;
    return min_value + t * (max_value - min_value);
}

static int random_int(int min_value, int max_value) {
    if (min_value >= max_value) {
        return min_value;
    }
    int range = max_value - min_value + 1;
    int value = rand() % range;
    return min_value + value;
}

static bool random_bool(double probability) {
    return random_double(0.0, 1.0) < probability;
}

typedef enum {
    PARAM_CLASS_INTEGER,
    PARAM_CLASS_RATIONAL,
    PARAM_CLASS_ENUM
} ParamClass;

typedef enum {
    PARAM_TARGET_TICKS,
    PARAM_TARGET_INITIAL_KOPPA,
    PARAM_TARGET_PSI_ALPHA,
    PARAM_TARGET_PSI_BETA,
    PARAM_TARGET_PSI_GAMMA,
    PARAM_TARGET_PSI_DELTA,
    PARAM_TARGET_PSI_EPSILON,
    PARAM_TARGET_KOPPA_ALPHA,
    PARAM_TARGET_KOPPA_BETA,
    PARAM_TARGET_KOPPA_GAMMA,
    PARAM_TARGET_KOPPA_DELTA,
    PARAM_TARGET_KOPPA_EPSILON,
    PARAM_TARGET_KOPPA_ZETA,
    PARAM_TARGET_KOPPA_ETA,
    PARAM_TARGET_KOPPA_THETA,
    PARAM_TARGET_KOPPA_IOTA,
    PARAM_TARGET_KOPPA_KAPPA,
    PARAM_TARGET_ENGINE_MODE,
    PARAM_TARGET_ENGINE_UPSILON,
    PARAM_TARGET_ENGINE_BETA,
    PARAM_TARGET_PSI_MODE,
    PARAM_TARGET_KOPPA_MODE,
    PARAM_TARGET_KOPPA_TRIGGER,
    PARAM_TARGET_PRIME_TARGET,
    PARAM_TARGET_MT10_BEHAVIOR,
    PARAM_TARGET_COUNT
} ParamTarget;

typedef struct {
    ParamClass class;
    ParamTarget target;
    double mutation_rate;
    double min_value;
    double max_value;
} ParamSpec;

static ParamSpec PARAM_SPECS[] = {
    {PARAM_CLASS_INTEGER, PARAM_TARGET_TICKS, 0.5, 10.0, 100.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_INITIAL_KOPPA, 0.5, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_PSI_ALPHA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_PSI_BETA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_PSI_GAMMA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_PSI_DELTA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_PSI_EPSILON, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_ALPHA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_BETA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_GAMMA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_DELTA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_EPSILON, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_ZETA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_ETA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_THETA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_IOTA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_RATIONAL, PARAM_TARGET_KOPPA_KAPPA, 0.3, -2.0, 2.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_ENGINE_MODE, 0.2, 0.0, 0.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_ENGINE_UPSILON, 0.2, 0.0, 0.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_ENGINE_BETA, 0.2, 0.0, 0.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_PSI_MODE, 0.2, 0.0, 0.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_KOPPA_MODE, 0.2, 0.0, 0.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_KOPPA_TRIGGER, 0.2, 0.0, 0.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_PRIME_TARGET, 0.2, 0.0, 0.0},
    {PARAM_CLASS_ENUM, PARAM_TARGET_MT10_BEHAVIOR, 0.2, 0.0, 0.0},
};

static EngineMode ENGINE_MODES[] = {
    ENGINE_MODE_ADD,
    ENGINE_MODE_MULTI,
    ENGINE_MODE_SLIDE,
    ENGINE_MODE_DELTA_ADD,
};

static EngineTrackMode ENGINE_TRACK_MODES[] = {
    ENGINE_TRACK_ADD,
    ENGINE_TRACK_MULTI,
    ENGINE_TRACK_SLIDE,
};

static PsiMode PSI_MODES[] = {
    PSI_MODE_STANDARD,
    PSI_MODE_ALTERNATE,
    PSI_MODE_EXPONENTIAL,
};

static KoppaMode KOPPA_MODES[] = {
    KOPPA_MODE_DUMP,
    KOPPA_MODE_POP,
    KOPPA_MODE_ACCUMULATE,
};

static KoppaTrigger KOPPA_TRIGGERS[] = {
    KOPPA_ON_ALL_MU,
    KOPPA_ON_PRIME_MU,
    KOPPA_ON_PERFECT_MU,
};

static PrimeTarget PRIME_TARGETS[] = {
    PRIME_ON_MEMORY,
    PRIME_ON_KOPPA,
    PRIME_ON_PSI,
};

static Mt10Behavior MT10_BEHAVIORS[] = {
    MT10_FORCED_PSI,
    MT10_FORCED_ENGINE,
    MT10_FORCED_KOPPA,
};

static EngineTrackMode evolution_track_mode_for_engine(EngineMode mode) {
    switch (mode) {
    case ENGINE_MODE_ADD:
        return ENGINE_TRACK_ADD;
    case ENGINE_MODE_MULTI:
        return ENGINE_TRACK_MULTI;
    case ENGINE_MODE_SLIDE:
        return ENGINE_TRACK_SLIDE;
    case ENGINE_MODE_DELTA_ADD:
    default:
        return ENGINE_TRACK_ADD;
    }
}

static void randomize_rational(Rational *value, double min_value, double max_value) {
    double num = random_double(min_value, max_value);
    double den = random_double(0.5, 2.0);
    rational_set_si(*value, (long)round(num * 1000.0), (unsigned long)round(den * 1000.0));
}

static void randomize_param(Config *config, ParamSpec *spec) {
    switch (spec->class) {
    case PARAM_CLASS_INTEGER: {
        if (spec->target == PARAM_TARGET_TICKS) {
            config->ticks = (unsigned)random_int((int)spec->min_value, (int)spec->max_value);
        }
    } break;
    case PARAM_CLASS_RATIONAL: {
        switch (spec->target) {
        case PARAM_TARGET_INITIAL_KOPPA:
            randomize_rational(&config->initial_koppa, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_PSI_ALPHA:
            randomize_rational(&config->psi_alpha, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_PSI_BETA:
            randomize_rational(&config->psi_beta, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_PSI_GAMMA:
            randomize_rational(&config->psi_gamma, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_PSI_DELTA:
            randomize_rational(&config->psi_delta, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_PSI_EPSILON:
            randomize_rational(&config->psi_epsilon, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_ALPHA:
            randomize_rational(&config->koppa_alpha, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_BETA:
            randomize_rational(&config->koppa_beta, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_GAMMA:
            randomize_rational(&config->koppa_gamma, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_DELTA:
            randomize_rational(&config->koppa_delta, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_EPSILON:
            randomize_rational(&config->koppa_epsilon, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_ZETA:
            randomize_rational(&config->koppa_zeta, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_ETA:
            randomize_rational(&config->koppa_eta, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_THETA:
            randomize_rational(&config->koppa_theta, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_IOTA:
            randomize_rational(&config->koppa_iota, spec->min_value, spec->max_value);
            break;
        case PARAM_TARGET_KOPPA_KAPPA:
            randomize_rational(&config->koppa_kappa, spec->min_value, spec->max_value);
            break;
        default:
            break;
        }
    } break;
    case PARAM_CLASS_ENUM: {
        switch (spec->target) {
        case PARAM_TARGET_ENGINE_MODE:
            config->engine_mode = ENGINE_MODES[rand() % ARRAY_COUNT(ENGINE_MODES)];
            config->engine_upsilon = evolution_track_mode_for_engine(config->engine_mode);
            config->engine_beta = evolution_track_mode_for_engine(config->engine_mode);
            break;
        case PARAM_TARGET_ENGINE_UPSILON:
            config->engine_upsilon = evolution_track_mode_for_engine(config->engine_mode);
            break;
        case PARAM_TARGET_ENGINE_BETA:
            config->engine_beta = evolution_track_mode_for_engine(config->engine_mode);
            break;
        case PARAM_TARGET_PSI_MODE:
            config->psi_mode = PSI_MODES[rand() % ARRAY_COUNT(PSI_MODES)];
            break;
        case PARAM_TARGET_KOPPA_MODE:
            config->koppa_mode = KOPPA_MODES[rand() % ARRAY_COUNT(KOPPA_MODES)];
            break;
        case PARAM_TARGET_KOPPA_TRIGGER:
            config->koppa_trigger = KOPPA_TRIGGERS[rand() % ARRAY_COUNT(KOPPA_TRIGGERS)];
            break;
        case PARAM_TARGET_PRIME_TARGET:
            config->prime_target = PRIME_TARGETS[rand() % ARRAY_COUNT(PRIME_TARGETS)];
            break;
        case PARAM_TARGET_MT10_BEHAVIOR:
            config->mt10_behavior = MT10_BEHAVIORS[rand() % ARRAY_COUNT(MT10_BEHAVIORS)];
            break;
        default:
            break;
        }
    } break;
    default:
        break;
    }
}

static void randomize_config(Config *config) {
    config_clear(config);
    config_init(config);

    for (size_t i = 0; i < ARRAY_COUNT(PARAM_SPECS); ++i) {
        randomize_param(config, &PARAM_SPECS[i]);
    }

    config->ticks = 30U;
    rational_set_si(config->initial_koppa, 1, 1);
    config->koppa_trigger = KOPPA_ON_ALL_MU;
    config->prime_target = PRIME_ON_MEMORY;
    config->mt10_behavior = MT10_FORCED_PSI;
}

static void mutate_integer(unsigned *value, double min_value, double max_value, double mutation_rate) {
    if (random_bool(mutation_rate)) {
        int delta = random_int(-5, 5);
        int new_value = (int)(*value) + delta;
        if (new_value < (int)min_value) new_value = (int)min_value;
        if (new_value > (int)max_value) new_value = (int)max_value;
        *value = (unsigned)new_value;
    }
}

static void mutate_rational(Rational *value, double min_value, double max_value, double mutation_rate) {
    if (random_bool(mutation_rate)) {
        double current_num = (double)rational_get_num(*value) / (double)rational_get_den(*value);
        double delta = random_double(-0.5, 0.5);
        double new_num = current_num + delta;
        if (new_num < min_value) new_num = min_value;
        if (new_num > max_value) new_num = max_value;
        rational_set_si(*value, (long)round(new_num * 1000.0), 1000UL);
    }
}

static void mutate_enum_int(int *value, int min_value, int max_value, double mutation_rate) {
    if (random_bool(mutation_rate)) {
        int choice = rand() % 3;
        if (choice == 0) {
            *value = random_int(min_value, max_value);
        } else if (choice == 1) {
            if (*value > min_value) {
                (*value)--;
            }
        } else {
            if (*value < max_value) {
                (*value)++;
            }
        }
    }
}

static void mutate_config(Config *config) {
    for (size_t i = 0; i < ARRAY_COUNT(PARAM_SPECS); ++i) {
        ParamSpec *spec = &PARAM_SPECS[i];
        switch (spec->class) {
        case PARAM_CLASS_INTEGER:
            if (spec->target == PARAM_TARGET_TICKS) {
                mutate_integer(&config->ticks, spec->min_value, spec->max_value, spec->mutation_rate);
            }
            break;
        case PARAM_CLASS_RATIONAL:
            switch (spec->target) {
            case PARAM_TARGET_INITIAL_KOPPA:
                mutate_rational(&config->initial_koppa, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_PSI_ALPHA:
                mutate_rational(&config->psi_alpha, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_PSI_BETA:
                mutate_rational(&config->psi_beta, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_PSI_GAMMA:
                mutate_rational(&config->psi_gamma, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_PSI_DELTA:
                mutate_rational(&config->psi_delta, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_PSI_EPSILON:
                mutate_rational(&config->psi_epsilon, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_ALPHA:
                mutate_rational(&config->koppa_alpha, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_BETA:
                mutate_rational(&config->koppa_beta, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_GAMMA:
                mutate_rational(&config->koppa_gamma, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_DELTA:
                mutate_rational(&config->koppa_delta, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_EPSILON:
                mutate_rational(&config->koppa_epsilon, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_ZETA:
                mutate_rational(&config->koppa_zeta, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_ETA:
                mutate_rational(&config->koppa_eta, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_THETA:
                mutate_rational(&config->koppa_theta, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_IOTA:
                mutate_rational(&config->koppa_iota, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            case PARAM_TARGET_KOPPA_KAPPA:
                mutate_rational(&config->koppa_kappa, spec->min_value, spec->max_value, spec->mutation_rate);
                break;
            default:
                break;
            }
            break;
        case PARAM_CLASS_ENUM:
            switch (spec->target) {
            case PARAM_TARGET_ENGINE_MODE: {
                int engine_index = 0;
                for (size_t j = 0; j < ARRAY_COUNT(ENGINE_MODES); ++j) {
                    if (ENGINE_MODES[j] == config->engine_mode) {
                        engine_index = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&engine_index, 0, (int)ARRAY_COUNT(ENGINE_MODES) - 1, spec->mutation_rate);
                config->engine_mode = ENGINE_MODES[engine_index];
                config->engine_upsilon = evolution_track_mode_for_engine(config->engine_mode);
                config->engine_beta = evolution_track_mode_for_engine(config->engine_mode);
            } break;
            case PARAM_TARGET_ENGINE_UPSILON: {
                int idx = 0;
                for (size_t j = 0; j < ARRAY_COUNT(ENGINE_TRACK_MODES); ++j) {
                    if (ENGINE_TRACK_MODES[j] == config->engine_upsilon) {
                        idx = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&idx, 0, (int)ARRAY_COUNT(ENGINE_TRACK_MODES) - 1, spec->mutation_rate);
                config->engine_upsilon = ENGINE_TRACK_MODES[idx];
            } break;
            case PARAM_TARGET_ENGINE_BETA: {
                int idx = 0;
                for (size_t j = 0; j < ARRAY_COUNT(ENGINE_TRACK_MODES); ++j) {
                    if (ENGINE_TRACK_MODES[j] == config->engine_beta) {
                        idx = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&idx, 0, (int)ARRAY_COUNT(ENGINE_TRACK_MODES) - 1, spec->mutation_rate);
                config->engine_beta = ENGINE_TRACK_MODES[idx];
            } break;
            case PARAM_TARGET_PSI_MODE: {
                int idx = 0;
                for (size_t j = 0; j < ARRAY_COUNT(PSI_MODES); ++j) {
                    if (PSI_MODES[j] == config->psi_mode) {
                        idx = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&idx, 0, (int)ARRAY_COUNT(PSI_MODES) - 1, spec->mutation_rate);
                config->psi_mode = PSI_MODES[idx];
            } break;
            case PARAM_TARGET_KOPPA_MODE: {
                int idx = 0;
                for (size_t j = 0; j < ARRAY_COUNT(KOPPA_MODES); ++j) {
                    if (KOPPA_MODES[j] == config->koppa_mode) {
                        idx = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&idx, 0, (int)ARRAY_COUNT(KOPPA_MODES) - 1, spec->mutation_rate);
                config->koppa_mode = KOPPA_MODES[idx];
            } break;
            case PARAM_TARGET_KOPPA_TRIGGER: {
                int idx = 0;
                for (size_t j = 0; j < ARRAY_COUNT(KOPPA_TRIGGERS); ++j) {
                    if (KOPPA_TRIGGERS[j] == config->koppa_trigger) {
                        idx = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&idx, 0, (int)ARRAY_COUNT(KOPPA_TRIGGERS) - 1, spec->mutation_rate);
                config->koppa_trigger = KOPPA_TRIGGERS[idx];
            } break;
            case PARAM_TARGET_PRIME_TARGET: {
                int idx = 0;
                for (size_t j = 0; j < ARRAY_COUNT(PRIME_TARGETS); ++j) {
                    if (PRIME_TARGETS[j] == config->prime_target) {
                        idx = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&idx, 0, (int)ARRAY_COUNT(PRIME_TARGETS) - 1, spec->mutation_rate);
                config->prime_target = PRIME_TARGETS[idx];
            } break;
            case PARAM_TARGET_MT10_BEHAVIOR: {
                int idx = 0;
                for (size_t j = 0; j < ARRAY_COUNT(MT10_BEHAVIORS); ++j) {
                    if (MT10_BEHAVIORS[j] == config->mt10_behavior) {
                        idx = (int)j;
                        break;
                    }
                }
                mutate_enum_int(&idx, 0, (int)ARRAY_COUNT(MT10_BEHAVIORS) - 1, spec->mutation_rate);
                config->mt10_behavior = MT10_BEHAVIORS[idx];
            } break;
            default:
                break;
            }
            break;
        default:
            break;
        }
    }
}

static double evaluate_candidate(Candidate *candidate, const char *summary_dir, bool use_summary_dir) {
    if (candidate->evaluated) {
        return candidate->score;
    }

    RunSummary summary;
    run_summary_init(&summary);

    SimulationResult result;
    simulation_result_init(&result);

    char summary_path[512] = {0};
    if (use_summary_dir) {
        snprintf(summary_path, sizeof(summary_path), "%s/candidate_%p_summary.json", summary_dir, (void *)candidate);
    }

    simulate(&candidate->config, &summary, &result);

    candidate->summary = summary;

    double score = 0.0;
    score += (double)summary.prime_events * 1.0;
    score += (double)summary.perfect_events * 1.0;
    score += (double)summary.psi_events * 0.5;
    score += (double)summary.koppa_events * 0.5;
    score -= (double)summary.memory_overflows * 10.0;

    candidate->score = score;
    candidate->evaluated = true;

    if (use_summary_dir) {
        FILE *file = fopen(summary_path, "w");
        if (file) {
            fprintf(file, "{\n");
            fprintf(file, "  \"score\": %.6f,\n", candidate->score);
            fprintf(file, "  \"prime_events\": %u,\n", summary.prime_events);
            fprintf(file, "  \"perfect_events\": %u,\n", summary.perfect_events);
            fprintf(file, "  \"psi_events\": %u,\n", summary.psi_events);
            fprintf(file, "  \"koppa_events\": %u,\n", summary.koppa_events);
            fprintf(file, "  \"memory_overflows\": %u,\n", summary.memory_overflows);
            fprintf(file, "  \"max_stack_depth\": %u,\n", summary.max_stack_depth);
            fprintf(file, "  \"ratio_mean\": %.6f,\n", summary.ratio_mean);
            fprintf(file, "  \"ratio_variance\": %.6f,\n", summary.ratio_variance);
            fprintf(file, "  \"stack_summary\": \"%s\"\n", summary.stack_summary);
            fprintf(file, "}\n");
            fclose(file);
        }
    }

    run_summary_clear(&summary);
    simulation_result_clear(&result);

    return candidate->score;
}

static int compare_candidates_by_score_desc(const void *a, const void *b) {
    const Candidate *ca = (const Candidate *)a;
    const Candidate *cb = (const Candidate *)b;
    if (ca->score < cb->score) return 1;
    if (ca->score > cb->score) return -1;
    return 0;
}

static void evolve_population(Candidate *population, unsigned population_size, unsigned elite_count, const char *summary_dir, bool use_summary_dir) {
    for (unsigned i = 0; i < population_size; ++i) {
        evaluate_candidate(&population[i], summary_dir, use_summary_dir);
    }

    qsort(population, population_size, sizeof(Candidate), compare_candidates_by_score_desc);

    for (unsigned i = elite_count; i < population_size; ++i) {
        candidate_clear(&population[i]);
    }

    for (unsigned i = elite_count; i < population_size; ++i) {
        unsigned parent_index = rand() % elite_count;
        config_clone(&population[i].config, &population[parent_index].config);
        run_summary_init(&population[i].summary);
        population[i].score = 0.0;
        population[i].evaluated = false;

        mutate_config(&population[i].config);
    }
}

static void print_candidate(const Candidate *candidate) {
    printf("Score: %.6f\n", candidate->score);
    printf("Ticks: %u\n", candidate->config.ticks);
    printf("Engine mode: %d\n", (int)candidate->config.engine_mode);
    printf("Engine upsilon: %d\n", (int)candidate->config.engine_upsilon);
    printf("Engine beta: %d\n", (int)candidate->config.engine_beta);
    printf("Psi mode: %d\n", (int)candidate->config.psi_mode);
    printf("Koppa mode: %d\n", (int)candidate->config.koppa_mode);
    printf("Koppa trigger: %d\n", (int)candidate->config.koppa_trigger);
    printf("Prime target: %d\n", (int)candidate->config.prime_target);
    printf("MT10 behavior: %d\n", (int)candidate->config.mt10_behavior);
}

static void parse_evolution_arguments(int argc, char **argv, EvolutionOptions *options) {
    options->generations = 10U;
    options->population = 8U;
    options->elite = 2U;
    options->seed = (unsigned)time(NULL);
    options->use_summary_dir = false;
    options->summary_dir[0] = '\0';

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--generations") == 0 && i + 1 < argc) {
            options->generations = (unsigned)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--population") == 0 && i + 1 < argc) {
            options->population = (unsigned)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--elite") == 0 && i + 1 < argc) {
            options->elite = (unsigned)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            options->seed = (unsigned)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--summary-dir") == 0 && i + 1 < argc) {
            strncpy(options->summary_dir, argv[++i], sizeof(options->summary_dir) - 1);
            options->summary_dir[sizeof(options->summary_dir) - 1] = '\0';
            options->use_summary_dir = true;
        }
    }

    if (options->population < 2U) {
        options->population = 2U;
    }
    if (options->elite == 0U || options->elite > options->population) {
        options->elite = 1U;
    }
}

int main(int argc, char **argv) {
    EvolutionOptions options;
    parse_evolution_arguments(argc, argv, &options);
    srand(options.seed);

    Candidate *population = malloc(sizeof(Candidate) * options.population);
    if (!population) {
        fprintf(stderr, "Failed to allocate population\n");
        return 1;
    }

    for (unsigned i = 0; i < options.population; ++i) {
        candidate_init(&population[i]);
        randomize_config(&population[i].config);
    }

    for (unsigned gen = 0; gen < options.generations; ++gen) {
        evolve_population(population, options.population, options.elite, options.summary_dir, options.use_summary_dir);
    }

    evolve_population(population, options.population, options.elite, options.summary_dir, options.use_summary_dir);

    Candidate best = population[0];
    print_candidate(&best);

    for (unsigned i = 0; i < options.population; ++i) {
        candidate_clear(&population[i]);
    }
    free(population);

    return 0;
}

