#ifndef TRTS_ENGINE_H
#define TRTS_ENGINE_H

#include "config.h"
#include "state.h"
#include <stdbool.h>

/* Perform a single engine update for the given microtick (1–11).
 * Returns true on success, false if a division by zero or invalid
 * operation occurred.  Updates upsilon and beta, deltas, triangles
 * and sign‑flip/polarity flags. */
bool engine_step(const Config *cfg, TRTS_State *st, int microtick);

#endif /* TRTS_ENGINE_H */

