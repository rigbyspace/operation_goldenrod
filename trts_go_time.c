// trts_go_time.c
/*
 * trts_go_time.c
 *
 * Standalone CLI simulator for the TRTS engine.
 *
 * - Runs the pure propagation engine only (no analysis, no pattern detection,
 *   no primality checks, no canonicalization, no GCD/factorization).
 * - Writes a simple CSV stream to stdout (or to a file if --output is given).
 * - Accepts minimal configuration via CLI: ticks and initial seeds for upsilon,
 *   beta and koppa as simple rational strings "N/D".
 *
 * Usage examples:
 *   ./trts_go_time --ticks 100 --ups 3/2 --beta 5/3 --koppa 1/1 > values.csv
 *   ./trts_go_time --ticks 50 --ups 1/1 --beta 1/1
 *
 * This program depends only on the engine, koppa, psi, state and rational
 * helpers in the codebase and intentionally omits any analysis logic.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <gmp.h>

#include "config.h"
#include "state.h"
#include "engine.h"
#include "koppa.h"
#include "psi.h"
#include "rational.h"

/* Helper: parse a rational string "N/D" into mpq_t using the project's
   rational_set_si primitive (which sets components without canonicalization).
   Returns true on success. */
static bool parse_seed(const char *text, mpq_t out) {
    if (!text || !out) return false;
    const char *slash = strchr(text, '/');
    if (!slash) return false;
    size_t nlen = (size_t)(slash - text);
    size_t dlen = strlen(slash + 1);
    if (nlen == 0 || dlen == 0 || nlen >= 128 || dlen >= 128) return false;

    char nb[128], db[128];
    memcpy(nb, text, nlen);
    nb[nlen] = '\0';
    memcpy(db, slash + 1, dlen);
    db[dlen] = '\0';

    char *endptr = NULL;
    long num = strtol(nb, &endptr, 10);
    if (!endptr || *endptr != '\0') return false;
    endptr = NULL;
    unsigned long den = strtoul(db, &endptr, 10);
    if (!endptr || *endptr != '\0' || den == 0UL) return false;

    rational_set_si(out, num, den);
    return true;
}

/* Print CSV header and a single values row (mirrors values.csv layout used elsewhere).
   We print numerator and denominator for upsilon, beta and koppa plus a few
   bookkeeping fields (koppa stack numerators/denominators and stack size). */
static void print_csv_header(FILE *f) {
    fprintf(f,
        "tick,mt,upsilon_num,upsilon_den,beta_num,beta_den,koppa_num,koppa_den,"
        "koppa_stack0_num,koppa_stack0_den,koppa_stack1_num,koppa_stack1_den,"
        "koppa_stack2_num,koppa_stack2_den,koppa_stack3_num,koppa_stack3_den,"
        "koppa_stack_size\n");
}

static void print_state_row(FILE *f, size_t tick, int microtick, const TRTS_State *s) {
    gmp_fprintf(f, "%zu,%d,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd,",
                tick, microtick,
                mpq_numref(s->upsilon), mpq_denref(s->upsilon),
                mpq_numref(s->beta), mpq_denref(s->beta),
                mpq_numref(s->koppa), mpq_denref(s->koppa));
    for (size_t i = 0; i < 4; ++i) {
        gmp_fprintf(f, "%Zd,%Zd,", mpq_numref(s->koppa_stack[i]), mpq_denref(s->koppa_stack[i]));
    }
    fprintf(f, "%zu\n", s->koppa_stack_size);
}

/* Minimal runtime configuration parsed from CLI */
typedef struct {
    size_t ticks;
    mpq_t upsilon_seed;
    mpq_t beta_seed;
    mpq_t koppa_seed;
    bool ups_set;
    bool beta_set;
    bool koppa_set;
    const char *output_path;
} RunConfig;

static void runconfig_init(RunConfig *rc) {
    rc->ticks = 30;
    rational_init(rc->upsilon_seed);
    rational_init(rc->beta_seed);
    rational_init(rc->koppa_seed);
    rc->ups_set = rc->beta_set = rc->koppa_set = false;
    rc->output_path = NULL;
}

static void runconfig_clear(RunConfig *rc) {
    rational_clear(rc->upsilon_seed);
    rational_clear(rc->beta_seed);
    rational_clear(rc->koppa_seed);
}

/* Simple CLI parser */
static bool parse_args(int argc, char **argv, RunConfig *rc) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--ticks") == 0 && i + 1 < argc) {
            rc->ticks = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--ups") == 0 && i + 1 < argc) {
            if (!parse_seed(argv[++i], rc->upsilon_seed)) return false;
            rc->ups_set = true;
        } else if (strcmp(argv[i], "--beta") == 0 && i + 1 < argc) {
            if (!parse_seed(argv[++i], rc->beta_seed)) return false;
            rc->beta_set = true;
        } else if (strcmp(argv[i], "--koppa") == 0 && i + 1 < argc) {
            if (!parse_seed(argv[++i], rc->koppa_seed)) return false;
            rc->koppa_set = true;
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            rc->output_path = argv[++i];
        } else if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            return false;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            return false;
        }
    }
    return true;
}

static void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [--ticks N] [--ups N/D] [--beta N/D] [--koppa N/D] [--output file]\n"
        "Runs the TRTS engine only and emits a CSV of raw state snapshots.\n"
        "Examples:\n"
        "  %s --ticks 100 --ups 3/2 --beta 5/3\n"
        "  %s --output run.csv --ups 1/1 --beta 1/1 --koppa 0/1\n",
        prog, prog, prog);
}

/* Main: compose Config from minimal defaults, run simulation loop, emit CSV rows.
   The simulation loop mirrors the microtick sequencing used in the rest of the
   codebase: 11 microticks per tick with phases E (1,4,7,10), M (2,5,8,11), R (others).
   No pattern or heavy evaluation is performed inside the loop. */
int main(int argc, char **argv) {
    RunConfig rc;
    runconfig_init(&rc);
    if (!parse_args(argc, argv, &rc)) {
        print_usage(argv[0]);
        runconfig_clear(&rc);
        return EXIT_FAILURE;
    }

    FILE *out = stdout;
    if (rc.output_path) {
        out = fopen(rc.output_path, "w");
        if (!out) {
            perror("opening output file");
            runconfig_clear(&rc);
            return EXIT_FAILURE;
        }
    }

    /* Build minimal Config and initial TRTS_State */
    Config cfg;
    config_init(&cfg);

    /* Keep defaults for other flags; only set initial rationals if provided */
    if (rc.ups_set) {
        rational_set(cfg.initial_upsilon, rc.upsilon_seed);
    }
    if (rc.beta_set) {
        rational_set(cfg.initial_beta, rc.beta_seed);
    }
    if (rc.koppa_set) {
        rational_set(cfg.initial_koppa, rc.koppa_seed);
    }

    /* Ensure multi-level koppa and other features default off for a clean run */
    cfg.multi_level_koppa = false;
    cfg.enable_modular_wrap = false;
    cfg.koppa_trigger = KOPPA_ON_ALL_MU;
    cfg.ticks = rc.ticks;

    TRTS_State state;
    state_init(&state);
    state_reset(&state, &cfg);

    /* Emit CSV header */
    print_csv_header(out);

    for (size_t tick = 1; tick <= cfg.ticks; ++tick) {
        for (int microtick = 1; microtick <= 11; ++microtick) {
            char phase;
            switch (microtick) {
                case 1: case 4: case 7: case 10: phase = 'E'; break;
                case 2: case 5: case 8: case 11: phase = 'M'; break;
                default: phase = 'R'; break;
            }

            switch (phase) {
                case 'E': {
                    /* E-phase: compute epsilon and apply engine_step */
                    rational_set(state.epsilon, state.upsilon);
                    (void)engine_step(&cfg, &state, microtick);
                    /* Intentionally do not run pattern detection or primality here */
                    state.rho_pending = false;
                    state.rho_latched = false;
                    break;
                }
                case 'M': {
                    /* Memory phase: decide psi firing using only lightweight flag checks
                       that do not run heavy evaluation. To remain pure and minimal we
                       only apply the psi transform when allowed by config (here: MSTEP
                       mode defaults) and stack rules. This is intentionally conservative. */
                    bool allow_stack = true;
                    bool psi_firable = (cfg.psi_mode == PSI_MODE_MSTEP);
                    if (psi_firable && allow_stack) {
                        (void)psi_transform(&cfg, &state);
                    }
                    koppa_accrue(&cfg, &state, state.psi_recent, true, microtick);
                    break;
                }
                case 'R': {
                    /* Reset phase: accrue koppa without psi */
                    koppa_accrue(&cfg, &state, false, false, microtick);
                    break;
                }
            }

            /* Emit raw snapshot row */
            print_state_row(out, tick, microtick, &state);
        }
    }

    /* Clean up */
    state_clear(&state);
    config_clear(&cfg);
    runconfig_clear(&rc);

    if (out != stdout) fclose(out);
    return EXIT_SUCCESS;
}

