// psi.h
#ifndef PSI_H
#define PSI_H

#include <stdbool.h>
#include <gmp.h>
#include "config.h"
#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Attempts to execute a psi transform on the given state.
 * Returns true if the transform succeeded, false otherwise.
 */
bool psi_transform(const Config *config, TRTS_State *state);

#ifdef __cplusplus
}
#endif

#endif // PSI_H
