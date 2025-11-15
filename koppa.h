/* koppa.h - TRTS Koppa (κ) Operations
 *
 * Manages the koppa register and 4-level stack.
 * Blueprint reference: koppa.h/koppa.c from Section 1.
 *
 * Koppa operations:
 * - DUMP: κ ← 0
 * - POP: κ ← ε
 * - ACCUMULATE: κ ← κ + ε
 *
 * Multi-level mode maintains a 4-entry FIFO stack and samples
 * specific stack entries at microticks 5 and 11.
 */

#ifndef TRTS_KOPPA_H
#define TRTS_KOPPA_H

#include "config.h"
#include "state.h"
#include <stdbool.h>

/* Update koppa according to configuration.
 *
 * Parameters:
 * - cfg: Configuration controlling trigger and mode
 * - st: State containing koppa, epsilon, and stack
 * - psi_fired: Whether psi just fired this microtick
 * - is_memory_step: True if current phase is M (memory)
 * - microtick: Current microtick number (1-11)
 *
 * Behavior:
 * - Determines if trigger condition is met based on koppa_trigger
 * - If triggered, pushes to stack (if multi-level) and performs operation
 * - Adds (υ + β) to κ after operation
 * - Updates koppa_sample from stack at specific microticks
 * - Maintains psi_recent flag for KOPPA_ON_MU_AFTER_PSI mode
 */
void koppa_accrue(const Config *cfg, TRTS_State *st, bool psi_fired, 
                  bool is_memory_step, int microtick);

#endif /* TRTS_KOPPA_H */
