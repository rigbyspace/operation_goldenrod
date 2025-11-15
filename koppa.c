#include "koppa.h"
#include "rational.h"

static void koppa_dump(TRTS_State *st) {
    rational_set_si(&st->koppa, 0, 1);
}

static void koppa_pop(TRTS_State *st) {
    rational_set(&st->koppa, &st->epsilon);
}

static void koppa_accumulate(TRTS_State *st) {
    Rational tmp;
    rational_init(&tmp);
    rational_add(&tmp, &st->koppa, &st->epsilon);
    rational_set(&st->koppa, &tmp);
    rational_clear(&tmp);
}

static void koppa_stack_push(TRTS_State *st, const Rational *val) {
    if (st->koppa_stack_size == 4) {
        for (size_t i = 1; i < 4; ++i) {
            rational_set(&st->koppa_stack[i - 1], &st->koppa_stack[i]);
        }
        rational_set(&st->koppa_stack[3], val);
    } else {
        rational_set(&st->koppa_stack[st->koppa_stack_size], val);
        st->koppa_stack_size += 1;
    }
}

static void koppa_update_sample(TRTS_State *st, int microtick, bool multi_level) {
    rational_set(&st->koppa_sample, &st->koppa);
    st->koppa_sample_index = -1;
    if (!multi_level) return;
    if (microtick == 11 && st->koppa_stack_size > 0) {
        rational_set(&st->koppa_sample, &st->koppa_stack[0]);
        st->koppa_sample_index = 0;
    } else if (microtick == 5 && st->koppa_stack_size > 2) {
        rational_set(&st->koppa_sample, &st->koppa_stack[2]);
        st->koppa_sample_index = 2;
    }
}

void koppa_accrue(const Config *cfg, TRTS_State *st, bool psi_fired, bool is_memory_step, int microtick) {
    bool trigger = false;
    switch (cfg->koppa_trigger) {
    case KOPPA_ON_PSI:         trigger = psi_fired; break;
    case KOPPA_ON_MU_AFTER_PSI: trigger = (is_memory_step && !psi_fired && st->psi_recent); break;
    case KOPPA_ON_ALL_MU:       trigger = is_memory_step; break;
    }
    if (!trigger) {
        if (!psi_fired && cfg->koppa_trigger != KOPPA_ON_ALL_MU) {
            st->psi_recent = st->psi_recent && (cfg->koppa_trigger == KOPPA_ON_MU_AFTER_PSI);
        }
        koppa_update_sample(st, microtick, cfg->multi_level_koppa);
        return;
    }
    if (cfg->multi_level_koppa) koppa_stack_push(st, &st->koppa);
    switch (cfg->koppa_mode) {
    case KOPPA_MODE_DUMP:  koppa_dump(st); break;
    case KOPPA_MODE_POP:   koppa_pop(st);  break;
    case KOPPA_MODE_ACCUMULATE: koppa_accumulate(st); break;
    }
    /* add upsilon + beta to koppa */
    Rational tmp;
    rational_init(&tmp);
    rational_add(&tmp, &st->upsilon, &st->beta);
    Rational sum;
    rational_init(&sum);
    rational_add(&sum, &st->koppa, &tmp);
    rational_set(&st->koppa, &sum);
    rational_clear(&tmp);
    rational_clear(&sum);
    if (cfg->koppa_trigger == KOPPA_ON_MU_AFTER_PSI) st->psi_recent = false;
    else st->psi_recent = psi_fired;
    koppa_update_sample(st, microtick, cfg->multi_level_koppa);
}

