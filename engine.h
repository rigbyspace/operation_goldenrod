// engine.h
#ifndef ENGINE_H
#define ENGINE_H

#include <stdbool.h>
#include "config.h"
#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Executes a single engine step for the given microtick.
 * Returns true if the step succeeded, false on error (e.g., division by zero).
 */
bool engine_step(const Config *config, TRTS_State *state, int microtick);

#ifdef __cplusplus
}
#endif

#endif // ENGINE_H
