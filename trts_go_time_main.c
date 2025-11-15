/* trts_go_time_main.c - Minimal TRTS Runner
 *
 * Stripped-down simulation for quick testing.
 * Pure propagation only, no heavy pattern detection.
 */

#include "config.h"
#include "state.h"
#include "engine.h"
#include "koppa.h"
#include "psi.h"
#include "rational.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static bool parse_rational(const char *text, mpq_t value) {
    if (!text) return false;
    const char *slash = strchr(text, '/');
    if (!slash) return false;
    
    char nb[128], db[128];
    size_t nl = (size_t)(slash - text);
    size_t dl = strlen(slash + 1);
    if (nl >= sizeof(nb) || dl >= sizeof(db)) return false;
    
    memcpy(nb, text, nl);
    nb[nl] = '\0';
    memcpy(db, slash + 1, dl);
    db[dl] = '\0';
    
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    
    if (mpz_set_str(n, nb, 10) != 0 || mpz_set_str(d, db, 10) != 0) {
        mpz_clear(n);
        mpz_clear(d);
        return false;
    }
    
    mpq_set_num(value, n);
    mpq_set_den(value, d);
    if (mpz_sgn(n) == 0) {
        mpz_set_ui(d, 0UL);
        mpq_set_den(value, d);
    }
    
    mpz_clear(n);
    mpz_clear(d);
    return true;
}

static void print_csv_header(FILE *f) {
    fprintf(f, "tick,mt,upsilon_num,upsilon_den,beta_num,beta_den,koppa_num,koppa_den\n");
}

static void print_state_row(FILE *f, size_t tick, int mt, const TRTS_State *s) {
    gmp_fprintf(f, "%zu,%d,%Zd,%Zd,%Zd,%Zd,%Zd,%Zd\n",
                tick, mt,
                mpq_numref(s->upsilon), mpq_denref(s->upsilon),
                mpq_numref(s->beta), mpq_denref(s->beta),
                mpq_numref(s->koppa), mpq_denref(s->koppa));
}

int main(int argc, char **argv) {
    Config cfg;
    config_init(&cfg);
    cfg.ticks = 30;
    
    FILE *out = stdout;
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--ticks") == 0 && i + 1 < argc) {
            cfg.ticks = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--ups") == 0 && i + 1 < argc) {
            if (!parse_rational(argv[++i], cfg.initial_upsilon)) {
                fprintf(stderr, "Invalid upsilon\n");
                config_clear(&cfg);
                return 1;
            }
        } else if (strcmp(argv[i], "--beta") == 0 && i + 1 < argc) {
            if (!parse_rational(argv[++i], cfg.initial_beta)) {
                fprintf(stderr, "Invalid beta\n");
                config_clear(&cfg);
                return 1;
            }
        } else if (strcmp(argv[i], "--koppa") == 0 && i + 1 < argc) {
            if (!parse_rational(argv[++i], cfg.initial_koppa)) {
                fprintf(stderr, "Invalid koppa\n");
                config_clear(&cfg);
                return 1;
            }
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            out = fopen(argv[++i], "w");
            if (!out) {
                perror("open output");
                config_clear(&cfg);
                return 1;
            }
        }
    }
    
    TRTS_State state;
    state_init(&state);
    state_reset(&state, &cfg);
    
    print_csv_header(out);
    
    for (size_t tick = 1; tick <= cfg.ticks; ++tick) {
        for (int mt = 1; mt <= 11; ++mt) {
            char phase = (mt == 1 || mt == 4 || mt == 7 || mt == 10) ? 'E' :
                        (mt == 2 || mt == 5 || mt == 8 || mt == 11) ? 'M' : 'R';
            
            switch (phase) {
                case 'E':
                    rational_set(&state.epsilon, &state.upsilon);
                    (void)engine_step(&cfg, &state, mt);
                    break;
                case 'M': {
                    bool psi_firable = (cfg.psi_mode == PSI_MODE_MSTEP);
                    if (psi_firable) {
                        (void)psi_transform(&cfg, &state);
                    }
                    koppa_accrue(&cfg, &state, state.psi_recent, true, mt);
                    break;
                }
                case 'R':
                    koppa_accrue(&cfg, &state, false, false, mt);
                    break;
            }
            
            print_state_row(out, tick, mt, &state);
        }
    }
    
    state_clear(&state);
    config_clear(&cfg);
    
    if (out != stdout) {
        fclose(out);
    }
    
    return 0;
}
