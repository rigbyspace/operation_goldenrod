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
#include "rational.h"

#define ARRAY_COUNT(arr) (sizeof(arr) / sizeof((arr)[0]))

typedef struct {
    Config config;
    RunSummary summary;
    double score;
    bool evaluated;
} Candidate;

typedef struct {
    size_t generations;
    size_t population;
    size_t elite;
    unsigned int seed;
    char strategy[32];
    char target_constant[32];
    bool save_output;
    char output_path[256];
} EvolutionOptions;

static EngineMode ENGINE_MODES[] = {
    ENGINE_MODE_ADD,
    ENGINE_MODE_MULTI,
    ENGINE_MODE_SLIDE,
    ENGINE_MODE_DELTA_ADD
};

static PsiMode PSI_MODES[] = {
    PSI_MODE_MSTEP,
    PSI_MODE_RHO_ONLY,
    PSI_MODE_MSTEP_RHO,
    PSI_MODE_INHIBIT_RHO
};

static KoppaMode KOPPA_MODES[] = {
    KOPPA_MODE_DUMP,
    KOPPA_MODE_POP,
    KOPPA_MODE_ACCUMULATE
};

static EngineTrackMode track_mode_for_engine(EngineMode mode) {
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

static void candidate_init(Candidate *candidate) {
    config_init(&candidate->config);
    run_summary_init(&candidate->summary);
    candidate->score = 0.0;
    candidate->evaluated = false;

    // Baseline tweaks for evolution search
    candidate->config.ticks = 30U;
    rational_set_si(candidate->config.initial_koppa, 1, 1);
    candidate->config.koppa_trigger = KOPPA_ON_ALL_MU;
    candidate->config.prime_target = PRIME_ON_MEMORY;
    candidate->config.mt10_behavior = MT10_FORCED_PSI;
}

static void candidate_clear(Candidate *candidate) {
    config_clear(&candidate->config);
    run_summary_clear(&candidate->summary);
}

static void config_clone(Config *dest, const Config *src) {
    dest->psi_mode = src->psi_mode;
    dest->koppa_mode = src->koppa_mode;
    dest->engine_mode = src->engine_mode;
    dest->engine_upsilon = src->engine_upsilon;
    dest->engine_beta = src->engine_beta;
    dest->dual_track_mode = src->dual_track_mode;
    dest->triple_psi_mode = src->triple_psi_mode;
    dest->multi_level_koppa = src->multi_level_koppa;
    dest->koppa_trigger = src->koppa_trigger;
    dest->prime_target = src->prime_target;
    dest->mt10_behavior = src->mt10_behavior;
    dest->ratio_trigger_mode = src->ratio_trigger_mode;
    dest->ticks = src->ticks;

    rational_set(dest->initial_upsilon, src->initial_upsilon);
    rational_set(dest->initial_beta, src->initial_beta);
    rational_set(dest->initial_koppa, src->initial_koppa);

    dest->enable_asymmetric_cascade = src->enable_asymmetric_cascade;
    dest->enable_conditional_triple_psi = src->enable_conditional_triple_psi;
    dest->enable_koppa_gated_engine = src->enable_koppa_gated_engine;
    dest->enable_delta_cross_propagation = src->enable_delta_cross_propagation;
    dest->enable_delta_koppa_offset = src->enable_delta_koppa_offset;
    dest->enable_ratio_threshold_psi = src->enable_ratio_threshold_psi;
    dest->enable_stack_depth_modes = src->enable_stack_depth_modes;
    dest->enable_epsilon_phi_triangle = src->enable_epsilon_phi_triangle;
    dest->enable_sign_flip = src->enable_sign_flip;
    dest->enable_modular_wrap = src->enable_modular_wrap;
    dest->enable_psi_strength_parameter = src->enable_psi_strength_parameter;
    dest->enable_ratio_snapshot_logging = src->enable_ratio_snapshot_logging;
    dest->enable_feedback_oscillator = src->enable_feedback_oscillator;
    dest->enable_fibonacci_gate = src->enable_fibonacci_gate;
    dest->sign_flip_mode = src->sign_flip_mode;
    dest->koppa_wrap_threshold = src->koppa_wrap_threshold;

    dest->enable_ratio_custom_range = src->enable_ratio_custom_range;
    rational_set(dest->ratio_custom_lower, src->ratio_custom_lower);
    rational_set(dest->ratio_custom_upper, src->ratio_custom_upper);

    dest->enable_twin_prime_trigger = src->enable_twin_prime_trigger;
    dest->enable_fibonacci_trigger = src->enable_fibonacci_trigger;
    dest->enable_perfect_power_trigger = src->enable_perfect_power_trigger;

    mpz_set(dest->modulus_bound, src->modulus_bound);
}

static void candidate_copy(Candidate *dest, const Candidate *src) {
    config_clone(&dest->config, &src->config);
    run_summary_copy(&dest->summary, &src->summary);
    dest->score = src->score;
    dest->evaluated = src->evaluated;
}

static long random_range(long min_value, long max_value) {
    long span = max_value - min_value + 1L;
    if (span <= 0L) {
        return min_value;
    }
    return min_value + (long)(rand() % (int)span);
}

static void mutate_seed(mpq_t value) {
    int choice = rand() % 4;
    if (choice == 0) {
        mpz_add_ui(mpq_numref(value), mpq_numref(value), 1UL);
    } else if (choice == 1) {
        mpz_sub_ui(mpq_numref(value), mpq_numref(value), 1UL);
    } else if (choice == 2) {
        unsigned long den = mpz_get_ui(mpq_denref(value));
        if (den > 1UL) {
            mpz_sub_ui(mpq_denref(value), mpq_denref(value), 1UL);
        }
    } else {
        mpz_add_ui(mpq_denref(value), mpq_denref(value), 1UL);
    }
}

static void randomize_config(Config *config) {
    config->engine_mode = ENGINE_MODES[rand() % ARRAY_COUNT(ENGINE_MODES)];
    config->engine_upsilon = track_mode_for_engine(config->engine_mode);
    config->engine_beta = track_mode_for_engine(config->engine_mode);
    config->psi_mode = PSI_MODES[rand() % ARRAY_COUNT(PSI_MODES)];
    config->koppa_mode = KOPPA_MODES[rand() % ARRAY_COUNT(KOPPA_MODES)];
    config->triple_psi_mode = (rand() % 2) != 0;
    config->multi_level_koppa = (rand() % 2) != 0;
    config->ticks = 25U + (size_t)(rand() % 10);

    long ups_num = random_range(1L, 8L);
    unsigned long ups_den = (unsigned long)random_range(1L, 8L);
    mpz_set_si(mpq_numref(config->initial_upsilon), ups_num);
    mpz_set_ui(mpq_denref(config->initial_upsilon), ups_den);

    long beta_num = random_range(1L, 8L);
    unsigned long beta_den = (unsigned long)random_range(1L, 8L);
    mpz_set_si(mpq_numref(config->initial_beta), beta_num);
    mpz_set_ui(mpq_denref(config->initial_beta), beta_den);

    // Keep koppa seed simple but non-zero to start
    rational_set_si(config->initial_koppa, 1, 1);
}

static void mutate_config(Config *config) {
    int mutations = 1 + rand() % 3;
    for (int i = 0; i < mutations; ++i) {
        int choice = rand() % 6;
        switch (choice) {
        case 0:
            config->engine_mode = ENGINE_MODES[rand() % ARRAY_COUNT(ENGINE_MODES)];
            config->engine_upsilon = track_mode_for_engine(config->engine_mode);
            config->engine_beta = track_mode_for_engine(config->engine_mode);
            break;
        case 1:
            config->psi_mode = PSI_MODES[rand() % ARRAY_COUNT(PSI_MODES)];
            break;
        case 2:
            config->koppa_mode = KOPPA_MODES[rand() % ARRAY_COUNT(KOPPA_MODES)];
            break;
        case 3:
            config->triple_psi_mode = !config->triple_psi_mode;
            break;
        case 4:
            mutate_seed(config->initial_upsilon);
            break;
        case 5:
        default:
            mutate_seed(config->initial_beta);
            break;
        }
    }
}

static double evaluate_candidate(Candidate *candidate, const EvolutionOptions *options) {
    if (candidate->evaluated) {
        return candidate->score;
    }

    RunSummary summary;
    run_summary_init(&summary);

    Config config;
    config_init(&config);
    config_clone(&config, &candidate->config);

    bool ok = simulate_and_analyze(&config, &summary);
    if (!ok) {
        candidate->score = -INFINITY;
        candidate->evaluated = true;
        run_summary_clear(&summary);
        config_clear(&config);
        return candidate->score;
    }

    double score = 0.0;

    double target_value = 0.0;
    bool has_target = analysis_constant_value(options->target_constant, &target_value);
    if (summary.ratio_defined && has_target) {
        double delta = fabs(summary.final_ratio_snapshot - target_value);
        score -= delta;
    }

    score += (double)summary.psi_events * 0.1;
    score += (double)summary.rho_events * 0.05;

    score -= summary.psi_spacing_stddev * 0.01;
    score -= summary.ratio_variance * 0.01;

    candidate->score = score;
    candidate->evaluated = true;
    run_summary_copy(&candidate->summary, &summary);

    run_summary_clear(&summary);
    config_clear(&config);
    return score;
}

static int compare_candidates(const void *a, const void *b) {
    const Candidate *ca = (const Candidate *)a;
    const Candidate *cb = (const Candidate *)b;
    if (cb->score > ca->score) return 1;
    if (cb->score < ca->score) return -1;
    return 0;
}

static void print_candidate_summary(const Candidate *candidate,
                                    size_t generation,
                                    size_t index) {
    printf("Generation %zu, Candidate %zu: score=%.6f, "
           "psi_events=%zu, rho_events=%zu, mu_zero=%zu, "
           "ratio=%.10f, stack_depth=%.3f\n",
           generation, index,
           candidate->score,
           candidate->summary.psi_events,
           candidate->summary.rho_events,
           candidate->summary.mu_zero_events,
           candidate->summary.final_ratio_snapshot,
           candidate->summary.average_stack_depth);
}

static void save_best_to_json(const Candidate *candidate,
                              const EvolutionOptions *options) {
    const char *path = options->output_path;
    if (!options->save_output || !path || path[0] == '\0') {
        return;
    }

    FILE *file = fopen(path, "w");
    if (!file) {
        perror("save_best_to_json");
        return;
    }

    const Config *config = &candidate->config;
    const RunSummary *summary = &candidate->summary;

    fprintf(file, "{\n");
    fprintf(file, "  \"score\": %.10f,\n", candidate->score);
    fprintf(file, "  \"engine_mode\": %d,\n", (int)config->engine_mode);
    fprintf(file, "  \"engine_upsilon\": %d,\n", (int)config->engine_upsilon);
    fprintf(file, "  \"engine_beta\": %d,\n", (int)config->engine_beta);
    fprintf(file, "  \"psi_mode\": %d,\n", (int)config->psi_mode);
    fprintf(file, "  \"koppa_mode\": %d,\n", (int)config->koppa_mode);
    fprintf(file, "  \"triple_psi_mode\": %s,\n", config->triple_psi_mode ? "true" : "false");
    fprintf(file, "  \"multi_level_koppa\": %s,\n", config->multi_level_koppa ? "true" : "false");
    fprintf(file, "  \"ticks\": %zu,\n", config->ticks);
    fprintf(file, "  \"final_ratio_snapshot\": %.10f,\n", summary->final_ratio_snapshot);
    fprintf(file, "  \"psi_events\": %zu,\n", summary->psi_events);
    fprintf(file, "  \"rho_events\": %zu,\n", summary->rho_events);
    fprintf(file, "  \"mu_zero_events\": %zu\n", summary->mu_zero_events);
    fprintf(file, "}\n");

    fclose(file);
}

static void parse_arguments(int argc, char **argv, EvolutionOptions *options) {
    options->generations = 10U;
    options->population = 8U;
    options->elite = 2U;
    options->seed = (unsigned int)time(NULL);
    snprintf(options->strategy, sizeof(options->strategy), "%s", "hill-climb");
    snprintf(options->target_constant, sizeof(options->target_constant), "%s", "rho");
    options->save_output = false;
    options->output_path[0] = '\0';

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--generations") == 0 && i + 1 < argc) {
            options->generations = (size_t)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--population") == 0 && i + 1 < argc) {
            options->population = (size_t)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--elite") == 0 && i + 1 < argc) {
            options->elite = (size_t)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            options->seed = (unsigned int)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--strategy") == 0 && i + 1 < argc) {
            snprintf(options->strategy, sizeof(options->strategy), "%s", argv[++i]);
        } else if (strcmp(argv[i], "--target") == 0 && i + 1 < argc) {
            snprintf(options->target_constant, sizeof(options->target_constant), "%s", argv[++i]);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            options->save_output = true;
            snprintf(options->output_path, sizeof(options->output_path), "%s", argv[++i]);
        }
    }

    if (options->elite == 0U || options->elite > options->population) {
        options->elite = 1U;
    }
}

int main(int argc, char **argv) {
    EvolutionOptions options;
    parse_arguments(argc, argv, &options);
    srand(options.seed);

    Candidate *population = malloc(sizeof(Candidate) * options.population);
    Candidate *next_population = malloc(sizeof(Candidate) * options.population);
    if (!population || !next_population) {
        fprintf(stderr, "Failed to allocate population buffers.\n");
        free(population);
        free(next_population);
        return 1;
    }

    for (size_t i = 0; i < options.population; ++i) {
        candidate_init(&population[i]);
        randomize_config(&population[i].config);
        population[i].evaluated = false;
    }

    for (size_t generation = 0; generation < options.generations; ++generation) {
        for (size_t i = 0; i < options.population; ++i) {
            evaluate_candidate(&population[i], &options);
        }

        qsort(population, options.population, sizeof(Candidate), compare_candidates);

        if (options.population > 0U) {
            print_candidate_summary(&population[0], generation, 0U);
        }

        for (size_t i = 0; i < options.population; ++i) {
            candidate_init(&next_population[i]);
        }

        size_t elite_count = options.elite;
        if (elite_count > options.population) {
            elite_count = options.population;
        }

        for (size_t i = 0; i < elite_count; ++i) {
            candidate_copy(&next_population[i], &population[i]);
        }

        for (size_t i = elite_count; i < options.population; ++i) {
            size_t parent_index = (size_t)(rand() % (int)elite_count);
            candidate_copy(&next_population[i], &population[parent_index]);
            mutate_config(&next_population[i].config);
            next_population[i].evaluated = false;
            next_population[i].score = 0.0;
        }

        for (size_t i = 0; i < options.population; ++i) {
            candidate_clear(&population[i]);
        }

        Candidate *temp = population;
        population = next_population;
        next_population = temp;
    }

    qsort(population, options.population, sizeof(Candidate), compare_candidates);
    if (options.population > 0U) {
        save_best_to_json(&population[0], &options);
        print_candidate_summary(&population[0], options.generations, 0U);
    }

    for (size_t i = 0; i < options.population; ++i) {
        candidate_clear(&population[i]);
    }
    free(population);
    free(next_population);
    return 0;
}

