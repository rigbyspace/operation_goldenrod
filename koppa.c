/* koppa.c - TRTS Koppa Implementation
 *
 * Implements koppa operations, stack management, and sampling logic.
 */

#include "koppa.h"
#include "rational.h"
#include "state.h"
#include "config.h"

/* Koppa operation: DUMP sets κ = 0 */
static void koppa_dump(TRTS_State *st) {
    rational_set_si(&st->koppa, 0, 1);
}

/* Koppa operation: POP sets κ = ε */
static void koppa_pop(TRTS_State *st) {
    rational_set(&st->koppa, &st->epsilon);
}

/* Koppa operation: ACCUMULATE sets κ = κ + ε */
static void koppa_accumulate(TRTS_State *st) {
    Rational tmp;
    rational_init(&tmp);
    rational_add(&tmp, &st->koppa, &st->epsilon);
    rational_set(&st->koppa, &tmp);
    rational_clear(&tmp);
}

/* Push value onto koppa stack (FIFO, max 4 entries) */
static void koppa_stack_push(TRTS_State *st, const Rational *val) {
    if (st->koppa_stack_size == 4) {
        /* Stack full: shift left and replace last entry */
        for (size_t i = 1; i < 4; ++i) {
            rational_set(&st->koppa_stack[i - 1], &st->koppa_stack[i]);
        }
        rational_set(&st->koppa_stack[3], val);
    } else {
        /* Stack not full: push to next available slot */
        rational_set(&st->koppa_stack[st->koppa_stack_size], val);
        st->koppa_stack_size++;
    }
}

/* Update koppa_sample from the stack based on microtick */
static void koppa_update_sample(TRTS_State *st, int microtick, bool multi_level) {
    if (!multi_level || st->koppa_stack_size == 0) {
        /* If not multi-level or stack is empty, sample is always kappa */
        rational_set(&st->koppa_sample, &st->koppa);
        return;
    }

    /* Sample stack at MT 5 and MT 11 */
    int index = -1;
    if (microtick == 5 && st->koppa_stack_size >= 1) {
        index = 0; /* Sample oldest value (index 0) */
    } else if (microtick == 11 && st->koppa_stack_size >= 3) {
        index = 2; /* Sample second to newest value (index 2) */
    }

    if (index != -1) {
        rational_set(&st->koppa_sample, &st->koppa_stack[index]);
        st->koppa_sample_index = index;
    } else {
        /* Otherwise, sample is the current koppa */
        rational_set(&st->koppa_sample, &st->koppa);
        st->koppa_sample_index = -1;
    }
}

/* Update koppa according to configuration. */
void koppa_accrue(const Config *cfg, TRTS_State *st, bool psi_fired, 
                  bool is_memory_step, int microtick) {
    
    /* 1. Determine if koppa operation should be triggered */
    bool trigger = false;
    switch (cfg->koppa_trigger) {
        case KOPPA_ON_PSI:
            trigger = psi_fired;
            break;
        case KOPPA_ON_MSTEP: // <-- FIX for warning: Added handler for KOPPA_ON_MSTEP
            trigger = is_memory_step;
            break;
        case KOPPA_ON_ALL_MU:
            trigger = true;
            break;
        case KOPPA_ON_MU_AFTER_PSI: // <-- FIX for error: Now declared in config.h
            trigger = st->psi_recent;
            break;
    }
    
    /* 2. Update psi_recent flag for KOPPA_ON_MU_AFTER_PSI mode */
    if (psi_fired) {
        st->psi_recent = true;
    } else if (cfg->koppa_trigger == KOPPA_ON_MU_AFTER_PSI && !is_memory_step) {
        st->psi_recent = true;
    } else {
        st->psi_recent = false;
    }

    if (!trigger) {
        /* If not triggered, just update the sample and exit. */
        koppa_update_sample(st, microtick, cfg->multi_level_koppa);
        return;
    }
    
    /* 3. Perform koppa operation based on mode */

    /* Push current koppa to stack if multi-level enabled */
    if (cfg->multi_level_koppa) {
        koppa_stack_push(st, &st->koppa);
    }
    
    switch (cfg->koppa_mode) {
        case KOPPA_MODE_DUMP:
            koppa_dump(st);
            break;
        case KOPPA_MODE_POP: // <-- FIX for error: Now declared in config.h
            koppa_pop(st);
            break;
        case KOPPA_MODE_ACCUMULATE:
            koppa_accumulate(st);
            break;
    }
    
    /* 4. Add (υ + β) to κ */
    Rational tmp, sum;
    rational_init(&tmp);
    rational_init(&sum);
    
    rational_add(&tmp, &st->upsilon, &st->beta);
    rational_add(&sum, &st->koppa, &tmp);
    rational_set(&st->koppa, &sum);
    
    /* 5. Cleanup and sample update */
    rational_clear(&tmp);
    rational_clear(&sum);
    
    koppa_update_sample(st, microtick, cfg->multi_level_koppa);
}
