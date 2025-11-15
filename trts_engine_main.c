// trts_engine_main.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "config_loader.h"
#include "simulate.h"

typedef struct {
    const Config *config;
} ObserverContext;

static const char *psi_mode_label(PsiMode mode) {
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

static void format_rational(mpq_srcptr value, char *buffer, size_t capacity) {
    if (!buffer || capacity == 0) {
        return;
    }

    char *num_str = mpz_get_str(NULL, 10, mpq_numref(value));
    char *den_str = mpz_get_str(NULL, 10, mpq_denref(value));
    if (!num_str || !den_str) {
        if (buffer && capacity > 0) {
            buffer[0] = '\0';
        }
    } else {
        snprintf(buffer, capacity, "%s/%s", num_str, den_str);
    }

    if (num_str) {
        free(num_str);
    }
    if (den_str) {
        free(den_str);
    }
}

static void append_token(char *buffer, size_t capacity, const char *token) {
    if (!buffer || capacity == 0 || !token || token[0] == '\0') {
        return;
    }

    size_t length = strlen(buffer);
    if (length >= capacity - 1) {
        return;
    }

    if (length > 0) {
        buffer[length++] = '|';
        if (length >= capacity - 1) {
            buffer[capacity - 1] = '\0';
            return;
        }
        buffer[length] = '\0';
    }

    strncat(buffer, token, capacity - strlen(buffer) - 1);
}

static void gui_observer(void *user_data, size_t tick, int microtick, char phase,
                         const TRTS_State *state, bool rho_event, bool psi_fired, bool mu_zero,
                         bool forced_emission) {
    const ObserverContext *context = (const ObserverContext *)user_data;

    char upsilon_buffer[256];
    char beta_buffer[256];
    char koppa_buffer[256];
    format_rational(state->upsilon, upsilon_buffer, sizeof(upsilon_buffer));
    format_rational(state->beta, beta_buffer, sizeof(beta_buffer));
    format_rational(state->koppa, koppa_buffer, sizeof(koppa_buffer));

    char psi_buffer[32];
    if (psi_fired) {
        if (state->psi_triple_recent) {
            snprintf(psi_buffer, sizeof(psi_buffer), "PSI_FIRE_TRIPLE");
        } else {
            snprintf(psi_buffer, sizeof(psi_buffer), "PSI_FIRE");
        }
    } else {
        snprintf(psi_buffer, sizeof(psi_buffer), "PSI_IDLE");
    }

    char rho_buffer[32];
    if (rho_event) {
        snprintf(rho_buffer, sizeof(rho_buffer), "RHO_EVENT");
    } else if (state->rho_pending || state->rho_latched) {
        snprintf(rho_buffer, sizeof(rho_buffer), "RHO_PENDING");
    } else {
        snprintf(rho_buffer, sizeof(rho_buffer), "RHO_IDLE");
    }

    char events_buffer[256];
    events_buffer[0] = '\0';
    char phase_token[2] = {phase, '\0'};
    append_token(events_buffer, sizeof(events_buffer), phase_token);
    if (mu_zero) {
        append_token(events_buffer, sizeof(events_buffer), "mu=0");
    }
    if (forced_emission) {
        append_token(events_buffer, sizeof(events_buffer), "forced");
    }
    if (state->ratio_triggered_recent) {
        append_token(events_buffer, sizeof(events_buffer), "ratio");
    }
    if (state->ratio_threshold_recent) {
        append_token(events_buffer, sizeof(events_buffer), "threshold");
    }
    if (state->dual_engine_last_step) {
        append_token(events_buffer, sizeof(events_buffer), "dual");
    }
    if (state->psi_strength_applied) {
        append_token(events_buffer, sizeof(events_buffer), "psi_strength");
    }
    if (state->koppa_sample_index >= 0) {
        char sample_token[32];
        snprintf(sample_token, sizeof(sample_token), "sample=%d", state->koppa_sample_index);
        append_token(events_buffer, sizeof(events_buffer), sample_token);
    }
    if (events_buffer[0] == '\0') {
        snprintf(events_buffer, sizeof(events_buffer), "-");
    }

    const char *psi_mode = context ? psi_mode_label(context->config->psi_mode) : "unknown";

    const char *color = psi_fired ? "#ff6f3c" : "#3a86ff";
    if (state->psi_triple_recent) {
        color = "#ff006e";
    } else if (state->ratio_triggered_recent) {
        color = "#00b894";
    } else if (state->ratio_threshold_recent) {
        color = "#8338ec";
    }

    printf("%zu;%d;%s;%s;%s;%s;%s;%zu;%s;%s;%s\n", tick, microtick, upsilon_buffer, beta_buffer,
           koppa_buffer, psi_buffer, rho_buffer, state->koppa_stack_size, events_buffer,
           psi_mode, color);
    fflush(stdout);
}

static void usage(const char *program) {
    fprintf(stderr, "Usage: %s --config <path>\n", program);
}

int main(int argc, char **argv) {
    const char *config_path = NULL;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--config") == 0) {
            if (i + 1 >= argc) {
                usage(argv[0]);
                return EXIT_FAILURE;
            }
            config_path = argv[++i];
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return EXIT_SUCCESS;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (!config_path) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    Config config;
    config_init(&config);

    char error_buffer[256];
    if (!config_load_from_file(&config, config_path, error_buffer, sizeof(error_buffer))) {
        fprintf(stderr, "Failed to load configuration: %s\n",
                (error_buffer[0] != '\0') ? error_buffer : "unknown error");
        config_clear(&config);
        return EXIT_FAILURE;
    }

    ObserverContext context = {&config};
    simulate_stream(&config, gui_observer, &context);

    config_clear(&config);
    return EXIT_SUCCESS;
}

