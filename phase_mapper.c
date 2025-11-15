// phase_mapper.c
#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulate.h"


#include "analysis_utils.h"
#include "config.h"

#define MAX_RESULTS 8192
#define ARRAY_COUNT(arr) (sizeof(arr) / sizeof((arr)[0]))

typedef struct {
    long numerator;
    unsigned long denominator;
} FractionSeed;

typedef struct {
    char engine[16];
    char psi[16];
    char koppa[16];
    char psi_type[16];
    char upsilon_seed[32];
    char beta_seed[32];
    char final_ratio[128];
    char closest_constant[32];
    char pattern[32];
    char classification[64];
    char stack_summary[128];
    double delta;
    size_t convergence_tick;
    double final_ratio_snapshot;
    double psi_spacing_mean;
    double psi_spacing_stddev;
    size_t psi_events;
    size_t rho_events;
    size_t mu_zero_events;
    double ratio_variance;
    double ratio_range;
    double ratio_stddev;
    double average_stack_depth;
} PhaseRecord;

typedef struct {
    bool scan_all;
    size_t ticks;
    size_t limit;
    bool verbose;
    bool write_output;
    char output_prefix[256];
    FractionSeed seeds[32];
    size_t seed_count;
} PhaseOptions;

static const char *engine_mode_name(EngineMode mode) {
    switch (mode) {
    case ENGINE_MODE_ADD:
        return "add";
    case ENGINE_MODE_MULTI:
        return "multi";
    case ENGINE_MODE_SLIDE:
        return "slide";
    case ENGINE_MODE_DELTA_ADD:
        return "delta";
    }
    return "unknown";
}

static const char *psi_mode_name(PsiMode mode) {
    switch (mode) {
    case PSI_MODE_MSTEP:
        return "mstep";
    case PSI_MODE_RHO_ONLY:
        return "rho_only";
    case PSI_MODE_MSTEP_RHO:
        return "mstep_rho";
    case PSI_MODE_INHIBIT_RHO:
        return "inhibit_rho";
    }
    return "unknown";
}

static const char *koppa_mode_name(KoppaMode mode) {
    switch (mode) {
    case KOPPA_MODE_DUMP:
        return "dump";
    case KOPPA_MODE_POP:
        return "pop";
    case KOPPA_MODE_ACCUMULATE:
        return "accumulate";
    }
    return "unknown";
}

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

static void config_to_strings(const Config *config, char *engine_out, size_t engine_size,
                              char *psi_out, size_t psi_size, char *koppa_out,
                              size_t koppa_size, char *psi_type_out, size_t psi_type_size) {
    snprintf(engine_out, engine_size, "%s", engine_mode_name(config->engine_mode));
    snprintf(psi_out, psi_size, "%s", psi_mode_name(config->psi_mode));
    snprintf(koppa_out, koppa_size, "%s", koppa_mode_name(config->koppa_mode));
    snprintf(psi_type_out, psi_type_size, "%s", analysis_psi_type_label(config));
}

static void apply_seed(mpq_t dest, FractionSeed seed) {
    mpz_set_si(mpq_numref(dest), seed.numerator);
    mpz_set_ui(mpq_denref(dest), seed.denominator);
}

static void options_init(PhaseOptions *options) {
    options->scan_all = false;
    options->ticks = 30U;
    options->limit = 0U;
    options->verbose = false;
    options->write_output = false;
    options->output_prefix[0] = '\0';
    options->seed_count = 0U;
}

static bool parse_fraction(const char *text, FractionSeed *seed) {
    if (!text || !seed) {
        return false;
    }
    const char *slash = strchr(text, '/');
    if (!slash) {
        return false;
    }
    char numerator_buffer[32];
    char denominator_buffer[32];
    size_t numerator_length = (size_t)(slash - text);
    size_t denominator_length = strlen(slash + 1);
    if (numerator_length >= sizeof(numerator_buffer) ||
        denominator_length >= sizeof(denominator_buffer)) {
        return false;
    }
    memcpy(numerator_buffer, text, numerator_length);
    numerator_buffer[numerator_length] = '\0';
    memcpy(denominator_buffer, slash + 1, denominator_length);
    denominator_buffer[denominator_length] = '\0';
    seed->numerator = strtol(numerator_buffer, NULL, 10);
    unsigned long denominator = strtoul(denominator_buffer, NULL, 10);
    if (denominator == 0UL) {
        return false;
    }
    seed->denominator = denominator;
    return true;
}

static void add_default_seeds(PhaseOptions *options) {
    static const FractionSeed defaults[] = {
        {1, 1}, {3, 2}, {5, 3}, {8, 5}, {7, 5}, {13, 8}
    };
    for (size_t i = 0; i < ARRAY_COUNT(defaults); ++i) {
        options->seeds[options->seed_count++] = defaults[i];
    }
}

static void parse_arguments(int argc, char **argv, PhaseOptions *options) {
    options_init(options);
    add_default_seeds(options);

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--scan-all") == 0) {
            options->scan_all = true;
        } else if (strcmp(argv[i], "--ticks") == 0 && i + 1 < argc) {
            options->ticks = (size_t)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--limit") == 0 && i + 1 < argc) {
            options->limit = (size_t)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--verbose") == 0) {
            options->verbose = true;
        } else if (strcmp(argv[i], "--output-phase-map") == 0 && i + 1 < argc) {
            options->write_output = true;
            snprintf(options->output_prefix, sizeof(options->output_prefix), "%s", argv[++i]);
        } else if (strcmp(argv[i], "--grid") == 0 && i + 1 < argc) {
            const char *grid_text = argv[++i];
            const char *delimiter = strchr(grid_text, ':');
            if (delimiter) {
                FractionSeed lower;
                FractionSeed upper;
                if (parse_fraction(grid_text, &lower) && parse_fraction(delimiter + 1, &upper)) {
                    options->seed_count = 0U;
                    for (long num = lower.numerator; num <= upper.numerator; ++num) {
                        for (unsigned long den = lower.denominator; den <= upper.denominator; ++den) {
                            if (options->seed_count >= ARRAY_COUNT(options->seeds)) {
                                break;
                            }
                            FractionSeed seed = {num, den};
                            options->seeds[options->seed_count++] = seed;
                        }
                    }
                }
            } else {
                options->seed_count = 0U;
                char buffer[256];
                snprintf(buffer, sizeof(buffer), "%s", grid_text);
                char *token = strtok(buffer, ",");
                while (token && options->seed_count < ARRAY_COUNT(options->seeds)) {
                    FractionSeed seed;
                    if (parse_fraction(token, &seed)) {
                        options->seeds[options->seed_count++] = seed;
                    }
                    token = strtok(NULL, ",");
                }
            }
            if (options->seed_count == 0U) {
                add_default_seeds(options);
            }
        }
    }
}

static void record_from_summary(const Config *config, const FractionSeed *ups_seed,
                                const FractionSeed *beta_seed, const RunSummary *summary,
                                PhaseRecord *record) {
    config_to_strings(config, record->engine, sizeof(record->engine), record->psi,
                      sizeof(record->psi), record->koppa, sizeof(record->koppa),
                      record->psi_type, sizeof(record->psi_type));

    snprintf(record->upsilon_seed, sizeof(record->upsilon_seed), "%ld/%lu", ups_seed->numerator,
             ups_seed->denominator);
    snprintf(record->beta_seed, sizeof(record->beta_seed), "%ld/%lu", beta_seed->numerator,
             beta_seed->denominator);

    snprintf(record->final_ratio, sizeof(record->final_ratio), "%s", summary->final_ratio_str);
    snprintf(record->closest_constant, sizeof(record->closest_constant), "%s",
             summary->closest_constant);
    snprintf(record->pattern, sizeof(record->pattern), "%s", summary->pattern);
    snprintf(record->classification, sizeof(record->classification), "%s", summary->classification);
    snprintf(record->stack_summary, sizeof(record->stack_summary), "%s", summary->stack_summary);

    record->delta = summary->closest_delta;
    record->convergence_tick = summary->convergence_tick;
    record->final_ratio_snapshot = summary->final_ratio_snapshot;
    record->psi_spacing_mean = summary->psi_spacing_mean;
    record->psi_spacing_stddev = summary->psi_spacing_stddev;
    record->psi_events = summary->psi_events;
    record->rho_events = summary->rho_events;
    record->mu_zero_events = summary->mu_zero_events;
    record->ratio_variance = summary->ratio_variance;
    record->ratio_range = summary->ratio_range;
    record->ratio_stddev = summary->ratio_stddev;
    record->average_stack_depth = summary->average_stack_depth;
}

static void write_csv(const PhaseRecord *records, size_t count, const char *path) {
    FILE *file = fopen(path, "w");
    if (!file) {
        return;
    }
    fprintf(file,
            "engine,psi,koppa,psi_type,u_seed,b_seed,final_ratio,closest_constant,delta,convergence_tick,pattern,classification,stack_summary,final_ratio_snapshot,psi_events,rho_events,mu_zero,psi_spacing_mean,psi_spacing_stddev,ratio_variance,ratio_range,ratio_stddev,average_stack_depth\n");
    for (size_t i = 0; i < count; ++i) {
        const PhaseRecord *record = &records[i];
        fprintf(file,
                "%s,%s,%s,%s,%s,%s,%s,%s,%.12g,%zu,%s,%s,%s,%.12g,%zu,%zu,%zu,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                record->engine, record->psi, record->koppa, record->psi_type, record->upsilon_seed,
                record->beta_seed, record->final_ratio, record->closest_constant, record->delta,
                record->convergence_tick, record->pattern, record->classification,
                record->stack_summary, record->final_ratio_snapshot, record->psi_events,
                record->rho_events, record->mu_zero_events, record->psi_spacing_mean,
                record->psi_spacing_stddev, record->ratio_variance, record->ratio_range,
                record->ratio_stddev, record->average_stack_depth);
    }
    fclose(file);
}

static void write_json(const PhaseRecord *records, size_t count, const char *path) {
    FILE *file = fopen(path, "w");
    if (!file) {
        return;
    }
    fprintf(file, "[\n");
    for (size_t i = 0; i < count; ++i) {
        const PhaseRecord *record = &records[i];
        fprintf(file,
                "  {\n"
                "    \"engine\": \"%s\",\n"
                "    \"psi\": \"%s\",\n"
                "    \"koppa\": \"%s\",\n"
                "    \"psi_type\": \"%s\",\n"
                "    \"upsilon_seed\": \"%s\",\n"
                "    \"beta_seed\": \"%s\",\n"
                "    \"final_ratio\": \"%s\",\n"
                "    \"closest_constant\": \"%s\",\n"
                "    \"delta\": %.12g,\n"
                "    \"convergence_tick\": %zu,\n"
                "    \"pattern\": \"%s\",\n"
                "    \"classification\": \"%s\",\n"
                "    \"stack_summary\": \"%s\",\n"
                "    \"final_ratio_snapshot\": %.12g,\n"
                "    \"psi_events\": %zu,\n"
                "    \"rho_events\": %zu,\n"
                "    \"mu_zero_events\": %zu,\n"
                "    \"psi_spacing_mean\": %.12g,\n"
                "    \"psi_spacing_stddev\": %.12g,\n"
                "    \"ratio_variance\": %.12g,\n"
                "    \"ratio_range\": %.12g,\n"
                "    \"ratio_stddev\": %.12g,\n"
                "    \"average_stack_depth\": %.12g\n"
                "  }%s\n",
                record->engine, record->psi, record->koppa, record->psi_type, record->upsilon_seed,
                record->beta_seed, record->final_ratio, record->closest_constant, record->delta,
                record->convergence_tick, record->pattern, record->classification,
                record->stack_summary, record->final_ratio_snapshot, record->psi_events,
                record->rho_events, record->mu_zero_events, record->psi_spacing_mean,
                record->psi_spacing_stddev, record->ratio_variance, record->ratio_range,
                record->ratio_stddev, record->average_stack_depth,
                (i + 1 < count) ? "," : "");
    }
    fprintf(file, "]\n");
    fclose(file);
}

static void print_record(const PhaseRecord *record) {
    printf("Engine=%s Psi=%s Koppa=%s Seeds=(%s,%s) Final=%s Constant=%s Δ=%.6g Pattern=%s Class=%s PSI=%s Stack=%s\n",
           record->engine, record->psi, record->koppa, record->upsilon_seed, record->beta_seed,
           record->final_ratio, record->closest_constant, record->delta, record->pattern,
           record->classification, record->psi_type, record->stack_summary);
}

int main(int argc, char **argv) {
    PhaseOptions options;
    parse_arguments(argc, argv, &options);

    if (options.seed_count == 0U) {
        fprintf(stderr, "No seeds available for phase mapping.\n");
        return 1;
    }

    PhaseRecord *records = malloc(sizeof(PhaseRecord) * MAX_RESULTS);
    if (!records) {
        fprintf(stderr, "Failed to allocate records.\n");
        return 1;
    }
    size_t record_count = 0U;

    Config config;
    config_init(&config);
    config.ticks = options.ticks;
    config.koppa_trigger = KOPPA_ON_ALL_MU;
    config.prime_target = PRIME_ON_MEMORY;
    config.mt10_behavior = MT10_FORCED_PSI;
    rational_set_si(config.initial_koppa, 1, 1);

    EngineMode engine_modes[] = {ENGINE_MODE_ADD, ENGINE_MODE_MULTI, ENGINE_MODE_SLIDE,
                                 ENGINE_MODE_DELTA_ADD};
    PsiMode psi_modes[] = {PSI_MODE_INHIBIT_RHO, PSI_MODE_MSTEP, PSI_MODE_RHO_ONLY,
                           PSI_MODE_MSTEP_RHO};
    KoppaMode koppa_modes[] = {KOPPA_MODE_DUMP, KOPPA_MODE_POP, KOPPA_MODE_ACCUMULATE};

    bool triple_modes[] = {false, true};

    for (size_t engine_index = 0; engine_index < ARRAY_COUNT(engine_modes); ++engine_index) {
        config.engine_mode = engine_modes[engine_index];
        config.engine_upsilon = track_mode_for_engine(config.engine_mode);
        config.engine_beta = track_mode_for_engine(config.engine_mode);

        for (size_t psi_index = 0; psi_index < ARRAY_COUNT(psi_modes); ++psi_index) {
            config.psi_mode = psi_modes[psi_index];

            for (size_t koppa_index = 0; koppa_index < ARRAY_COUNT(koppa_modes); ++koppa_index) {
                config.koppa_mode = koppa_modes[koppa_index];

                for (size_t triple_index = 0; triple_index < ARRAY_COUNT(triple_modes); ++triple_index) {
                    config.triple_psi_mode = triple_modes[triple_index];

                    for (size_t u_index = 0; u_index < options.seed_count; ++u_index) {
                        FractionSeed ups_seed = options.seeds[u_index];
                        apply_seed(config.initial_upsilon, ups_seed);

                        for (size_t b_index = 0; b_index < options.seed_count; ++b_index) {
                            FractionSeed beta_seed = options.seeds[b_index];
                            apply_seed(config.initial_beta, beta_seed);

                            RunSummary summary;
                            run_summary_init(&summary);
                            bool ok = simulate_and_analyze(&config, &summary);
                            if (!ok) {
                                run_summary_clear(&summary);
                                continue;
                            }

                            if (record_count < MAX_RESULTS) {
                                record_from_summary(&config, &ups_seed, &beta_seed, &summary,
                                                    &records[record_count]);
                                if (options.verbose) {
                                    print_record(&records[record_count]);
                                }
                                ++record_count;
                            }

                            if (options.limit > 0U && record_count >= options.limit) {
                                run_summary_clear(&summary);
                                goto phase_done;
                            }

                            run_summary_clear(&summary);
                        }
                    }
                }
            }
        }
    }

phase_done:
    if (options.write_output && record_count > 0U) {
        char csv_path[512];
        char json_path[512];
        snprintf(csv_path, sizeof(csv_path), "%s.csv", options.output_prefix);
        snprintf(json_path, sizeof(json_path), "%s.json", options.output_prefix);
        write_csv(records, record_count, csv_path);
        write_json(records, record_count, json_path);
    }

    if (!options.verbose) {
        for (size_t i = 0; i < record_count; ++i) {
            print_record(&records[i]);
        }
    }

    free(records);
    config_clear(&config);
    return 0;
}

