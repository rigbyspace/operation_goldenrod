// koppa.c
/*
==============================================
     TRTS SYSTEM CREED – RATIONAL ONLY
==============================================
*/

#include "koppa.h"

#include <gmp.h>

#include "rational.h"

static void koppa_dump(TRTS_State *state) {
    rational_set_si(state->koppa, 0, 1);
}

static void koppa_pop(TRTS_State *state) {
    rational_set(state->koppa, state->epsilon);
}

static void koppa_accumulate(TRTS_State *state) {
    rational_add(state->koppa, state->koppa, state->epsilon);
}

static void koppa_stack_push(TRTS_State *state, mpq_srcptr value) {
    if (state->koppa_stack_size == 4) {
        for (size_t i = 1; i < 4; ++i) {
            rational_set(state->koppa_stack[i - 1], state->koppa_stack[i]);
        }
        rational_set(state->koppa_stack[3], value);
    } else {
        rational_set(state->koppa_stack[state->koppa_stack_size], value);
        state->koppa_stack_size += 1;
    }
}

static void koppa_update_sample(TRTS_State *state, int microtick, bool multi_level_active) {
    state->koppa_sample_index = -1;
    rational_set(state->koppa_sample, state->koppa);

    if (!multi_level_active) {
        return;
    }

    if (microtick == 11 && state->koppa_stack_size > 0) {
        rational_set(state->koppa_sample, state->koppa_stack[0]);
        state->koppa_sample_index = 0;
    } else if (microtick == 5 && state->koppa_stack_size > 2) {
        rational_set(state->koppa_sample, state->koppa_stack[2]);
        state->koppa_sample_index = 2;
    }
}

void koppa_accrue(const Config *config, TRTS_State *state, bool psi_fired, bool is_memory_step,
                  int microtick) {
    bool trigger = false;

    switch (config->koppa_trigger) {
    case KOPPA_ON_PSI:
        trigger = psi_fired;
        break;
    case KOPPA_ON_MU_AFTER_PSI:
        trigger = is_memory_step && !psi_fired && state->psi_recent;
        break;
    case KOPPA_ON_ALL_MU:
        trigger = is_memory_step;
        break;
    }

    if (!trigger) {
        if (!psi_fired && config->koppa_trigger != KOPPA_ON_ALL_MU) {
            state->psi_recent = state->psi_recent && (config->koppa_trigger == KOPPA_ON_MU_AFTER_PSI);
        }
        koppa_update_sample(state, microtick, config->multi_level_koppa);
        return;
    }

    if (config->multi_level_koppa) {
        koppa_stack_push(state, state->koppa);
    }

    switch (config->koppa_mode) {
    case KOPPA_MODE_DUMP:
        koppa_dump(state);
        break;
    case KOPPA_MODE_POP:
        koppa_pop(state);
        break;
    case KOPPA_MODE_ACCUMULATE:
        koppa_accumulate(state);
        break;
    }

    mpq_t addition;
    rational_init(addition);
    rational_add(addition, state->upsilon, state->beta);
    rational_add(state->koppa, state->koppa, addition);
    rational_clear(addition);

    if (config->koppa_trigger == KOPPA_ON_MU_AFTER_PSI) {
        state->psi_recent = false;
    } else {
        state->psi_recent = psi_fired;
    }

    koppa_update_sample(state, microtick, config->multi_level_koppa);
}

