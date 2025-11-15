/* simulate.h - TRTS Simulation Loop
 *
 * Orchestrates the complete TRTS simulation with microtick sequencing,
 * pattern detection, ratio triggers, and event emission.
 * Blueprint reference: simulate.h/simulate.c from Section 1.
 *
 * Microtick phases (11 per tick):
 * - E (Epsilon): 1, 4, 7, 10 - compute ε and run engine step
 * - M (Memory): 2, 5, 8, 11 - check patterns, fire ψ, accrue κ
 * - R (Reset): 3, 6, 9 - accrue κ without ψ
 *
 * Emission events occur at microticks 1, 4, 7, 10 (E-phase).
 * Koppa acts as phantom microtick 12 continuation.
 */

#ifndef TRTS_SIMULATE_H
#define TRTS_SIMULATE_H

#include "config.h"
#include "state.h"
#include <stddef.h>
#include <stdbool.h>

/* Observer callback for streaming simulation data
 *
 * Called once per microtick with current state and event flags.
 * This callback is invoked DURING propagation, so the state contains
 * raw, non-canonicalized rationals.
 *
 * Parameters:
 * - user_data: User-provided context pointer
 * - tick: Current tick number (1-based)
 * - microtick: Current microtick (1-11)
 * - phase: Phase character ('E', 'M', or 'R')
 * - state: Current TRTS state (read-only)
 * - rho_event: True if ρ was triggered this microtick
 * - psi_fired: True if ψ transform was applied
 * - mu_zero: True if β numerator is zero
 * - forced_emission: True if this is microtick 10
 */
typedef void (*SimulateObserver)(void *user_data, size_t tick, int microtick, 
                                  char phase, const TRTS_State *state, 
                                  bool rho_event, bool psi_fired, 
                                  bool mu_zero, bool forced_emission);

/* Run simulation and write events.csv and values.csv
 *
 * This function executes a complete TRTS simulation with file output.
 * CSV files contain raw numerator/denominator values (no canonicalization).
 *
 * Files created:
 * - events.csv: Event flags per microtick
 * - values.csv: Raw rational values per microtick
 */
void simulate(const Config *config);

/* Run simulation with observer callback (no file output)
 *
 * This function executes a complete TRTS simulation and invokes
 * the observer callback at each microtick. Used for in-memory
 * analysis and interactive displays.
 *
 * Parameters:
 * - config: Simulation configuration
 * - observer: Callback function (may be NULL)
 * - user_data: Context pointer passed to observer
 */
void simulate_stream(const Config *config, SimulateObserver observer, 
                     void *user_data);

#endif /* TRTS_SIMULATE_H */
