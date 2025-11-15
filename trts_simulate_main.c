/* trts_simulate_main.c - Simple TRTS Simulator
 *
 * Runs TRTS simulation and writes events.csv and values.csv.
 */

#include "config.h"
#include "simulate.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(const char *prog) {
    fprintf(stderr, 
        "Usage: %s [options]\n"
        "Options:\n"
        "  --ticks N           Number of ticks to simulate (default: 10)\n"
        "  --ups N/D           Initial upsilon seed (default: 1/1)\n"
        "  --beta N/D          Initial beta seed (default: 1/1)\n"
        "  --koppa N/D         Initial koppa seed (default: 0/0)\n"
        "  --engine-mode N     Engine mode 0-3 (default: 0=ADD)\n"
        "  --psi-mode N        Psi mode 0-3 (default: 0=MSTEP)\n"
        "  --triple-psi        Enable 3-way psi transform\n"
        "  --multi-level       Enable 4-level koppa stack\n"
        "  -h, --help          Show this help\n\n"
        "Outputs: events.csv and values.csv\n",
        prog);
}

/* Parse rational string "N/D" */
static bool parse_rational(const char *text, mpq_t value) {
    if (!text || !value) {
        return false;
    }
    
    const char *slash = strchr(text, '/');
    if (!slash) {
        return false;
    }
    
    char num_buf[256], den_buf[256];
    size_t num_len = (size_t)(slash - text);
    size_t den_len = strlen(slash + 1);
    
    if (num_len >= sizeof(num_buf) || den_len >= sizeof(den_buf)) {
        return false;
    }
    
    memcpy(num_buf, text, num_len);
    num_buf[num_len] = '\0';
    memcpy(den_buf, slash + 1, den_len);
    den_buf[den_len] = '\0';
    
    mpz_t num, den;
    mpz_init(num);
    mpz_init(den);
    
    if (mpz_set_str(num, num_buf, 10) != 0 || mpz_set_str(den, den_buf, 10) != 0) {
        mpz_clear(num);
        mpz_clear(den);
        return false;
    }
    
    mpq_set_num(value, num);
    mpq_set_den(value, den);
    
    /* Enforce 0/0 for zero numerator */
    if (mpz_sgn(num) == 0) {
        mpz_set_ui(den, 0UL);
        mpq_set_den(value, den);
    }
    
    mpz_clear(num);
    mpz_clear(den);
    return true;
}

int main(int argc, char **argv) {
    Config config;
    config_init(&config);
    
    /* Parse arguments */
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--ticks") == 0 && i + 1 < argc) {
            config.ticks = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--ups") == 0 && i + 1 < argc) {
            if (!parse_rational(argv[++i], config.initial_upsilon)) {
                fprintf(stderr, "Invalid upsilon seed\n");
                config_clear(&config);
                return 1;
            }
        } else if (strcmp(argv[i], "--beta") == 0 && i + 1 < argc) {
            if (!parse_rational(argv[++i], config.initial_beta)) {
                fprintf(stderr, "Invalid beta seed\n");
                config_clear(&config);
                return 1;
            }
        } else if (strcmp(argv[i], "--koppa") == 0 && i + 1 < argc) {
            if (!parse_rational(argv[++i], config.initial_koppa)) {
                fprintf(stderr, "Invalid koppa seed\n");
                config_clear(&config);
                return 1;
            }
        } else if (strcmp(argv[i], "--engine-mode") == 0 && i + 1 < argc) {
            int mode = atoi(argv[++i]);
            if (mode >= 0 && mode <= 3) {
                config.engine_mode = (EngineMode)mode;
            }
        } else if (strcmp(argv[i], "--psi-mode") == 0 && i + 1 < argc) {
            int mode = atoi(argv[++i]);
            if (mode >= 0 && mode <= 3) {
                config.psi_mode = (PsiMode)mode;
            }
        } else if (strcmp(argv[i], "--triple-psi") == 0) {
            config.triple_psi_mode = true;
        } else if (strcmp(argv[i], "--multi-level") == 0) {
            config.multi_level_koppa = true;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            config_clear(&config);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            config_clear(&config);
            return 1;
        }
    }
    
    printf("TRTS Simulation\n");
    printf("===============\n");
    printf("Ticks: %zu\n", config.ticks);
    printf("Engine mode: %d\n", (int)config.engine_mode);
    printf("Psi mode: %d\n", (int)config.psi_mode);
    printf("Triple psi: %s\n", config.triple_psi_mode ? "yes" : "no");
    printf("Multi-level koppa: %s\n", config.multi_level_koppa ? "yes" : "no");
    printf("\nRunning simulation...\n");
    
    simulate(&config);
    
    printf("Complete. Output written to events.csv and values.csv\n");
    
    config_clear(&config);
    return 0;
}
