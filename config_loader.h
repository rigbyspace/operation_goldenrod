/* config_loader.h - TRTS Configuration Loader
 *
 * Parses JSON configuration files and populates Config structures.
 * Blueprint reference: config_loader.h/config_loader.c from Section 1.
 *
 * Supports loading:
 * - Mode enumerations (engine, psi, koppa)
 * - Boolean feature flags
 * - Integer parameters (ticks, thresholds)
 * - Rational seeds (as "N/D" strings)
 * - Custom ratio ranges
 * - Modulus bounds
 */

#ifndef TRTS_CONFIG_LOADER_H
#define TRTS_CONFIG_LOADER_H

#include "config.h"
#include <stdbool.h>
#include <stddef.h>

/* Load configuration from JSON file.
 *
 * Parameters:
 * - config: Config structure to populate
 * - path: Path to JSON file
 * - error_buffer: Buffer for error messages (may be NULL)
 * - error_capacity: Size of error buffer
 *
 * Returns: true on success, false on error
 *
 * JSON format example:
 * {
 *   "tick_count": 100,
 *   "engine_mode": 0,
 *   "psi_mode": 1,
 *   "upsilon_seed": "3/2",
 *   "beta_seed": "5/3",
 *   "triple_psi": true
 * }
 */
bool config_load_from_file(Config *config, const char *path, 
                           char *error_buffer, size_t error_capacity);

#endif /* TRTS_CONFIG_LOADER_H */
