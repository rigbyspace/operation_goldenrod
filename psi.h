#ifndef TRTS_PSI_H
#define TRTS_PSI_H

#include "config.h"
#include "state.h"
#include <stdbool.h>

/* Execute the psi transform if conditions allow.  Returns true if a transform
 * was applied.  Updates state->upsilon/beta/koppa and clears rho_pending on
 * the first successful fire. */
bool psi_transform(const Config *cfg, TRTS_State *st);

#endif /* TRTS_PSI_H */

