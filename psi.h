/* psi.h - TRTS Psi (ψ) Transform
 *
 * Implements the ψ transformation that exchanges and inverts rationals.
 * Blueprint reference: psi.h/psi.c from Section 1.
 *
 * Standard 2-way ψ: (υ,β) → (β/υ, υ/β)
 * Triple 3-way ψ: (υ,β,κ) → (β/κ, κ/υ, κ/β)
 *
 * All divisions are checked for zero denominators. No canonicalization.
 */

#ifndef TRTS_PSI_H
#define TRTS_PSI_H

#include "config.h"
#include "state.h"
#include <stdbool.h>

/* Execute psi transform if conditions allow.
 *
 * Behavior:
 * - Checks psi_mode and rho_pending to determine if firing is allowed
 * - Executes standard or triple psi based on configuration
 * - May fire multiple times if psi_strength_parameter is enabled
 * - Updates state flags (psi_recent, psi_triple_recent, psi_strength_applied)
 * - Clears rho_pending on first successful fire
 *
 * Returns: true if at least one transform was applied, false otherwise
 */
bool psi_transform(const Config *cfg, TRTS_State *st);

#endif /* TRTS_PSI_H */
