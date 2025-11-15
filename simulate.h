// simulate.h
#ifndef SIMULATE_H
#define SIMULATE_H

#include <stddef.h>
#include <stdbool.h>
#include "config.h"
#include "state.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Observer callback for streaming simulation results.
 * Called once per microtick with current state and event flags.
 */
typedef void (*SimulateObserver)(void *user_data, size_t tick, int microtick, char phase,
                                  const TRTS_State *state, bool rho_event, bool psi_fired,
                                  bool mu_zero, bool forced_emission);

/**
 * Runs a complete simulation and writes results to events.csv and values.csv.
 */
void simulate(const Config *config);

/**
 * Runs a simulation with a custom observer callback instead of file output.
 */
void simulate_stream(const Config *config, SimulateObserver observer, void *user_data);

#ifdef __cplusplus
}
#endif

#endif // SIMULATE_H

