/* engine.h - TRTS Engine Step
 *
 * Core propagation logic for updating υ and β.
 * Blueprint reference: engine.h/engine.c from Section 1.
 *
 * The engine step:
 * 1. Stores previous values of υ and β
 * 2. Determines track modes (ADD/MULTI/SLIDE) for υ and β
 * 3. Applies optional modulations (asymmetric cascade, stack depth, koppa gate)
 * 4. Computes new υ and β using track formulas
 * 5. Applies cross-propagation, sign-flip, and modular wrap
 * 6. Updates triangles, deltas, and state flags
 *
 * Returns false if a division by zero occurs, true otherwise.
 */

#ifndef TRTS_ENGINE_H
#define TRTS_ENGINE_H

#include "config.h"
#include "state.h"
#include <stdbool.h>

/* Perform a single engine update for the given microtick (1-11).
 *
 * Parameters:
 * - config: Configuration controlling engine behavior
 * - state: State containing υ, β, κ and other registers
 * - microtick: Current microtick number (affects asymmetric modes)
 *
 * Returns: true on success, false if division by zero or invalid operation
 *
 * Side effects:
 * - Updates state->upsilon and state->beta
 * - Updates state->previous_upsilon, state->previous_beta
 * - Updates state->delta_upsilon, state->delta_beta
 * - Updates triangle ratios if enabled
 * - Sets state->dual_engine_last_step flag
 * - Applies modular wrap if enabled
 */
bool engine_step(const Config *config, TRTS_State *state, int microtick);

#endif /* TRTS_ENGINE_H */
