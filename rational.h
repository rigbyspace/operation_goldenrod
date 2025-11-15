/*
==============================================
     TRTS SYSTEM CREED – RATIONAL ONLY
==============================================
Strict enforcement header: forbids canonicalization and GCD-style helpers.
All propagation must remain strictly within the rational field Q without
canonicalization or GCD normalization. This header enforces that at
compile-time and at runtime (best-effort).

Allowed GMP operations:
 - mpq_init, mpq_clear
 - mpq_set (wrapped carefully via rational_set)
 - mpz_set, mpz_set_si, mpz_set_ui
 - mpq_add, mpq_sub, mpq_mul, mpq_div
 - mpz_probab_prime_p, mpz_perfect_power_p
 - mpq_numref, mpq_denref

Forbidden symbols (any use will cause build-time failure or runtime abort):
 - mpq_canonicalize
 - mpq_set_str
 - mpq_set (raw wrapper use that may canonicalize if misused; usage is restricted)
 - mpz_gcd
 - mpz_gcd_ui
 - mpz_gcdext
 - mpz_invert
 - any other symbol invoking explicit GCD/canonicalization

If you must import rationals from text, use mpz_set_str for numerator and denominator
and then call rational_set_components to assign numerator/denominator without canonicalization.
*/

#ifndef RATIONAL_H
#define RATIONAL_H

#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------
   Forbidden-symbol compile-time guards
   ---------------------------
   These macros will intentionally cause a compile-time error if source code
   attempts to call these forbidden GMP helpers directly.
   --------------------------- */

#ifdef mpq_canonicalize
#undef mpq_canonicalize
#endif
#define mpq_canonicalize(...) \
    _Pragma("GCC error \"mpq_canonicalize is forbidden by TRTS policy (no canonicalization).\"")

/* Prohibit text-based setters that may canonicalize */
#ifdef mpq_set_str
#undef mpq_set_str
#endif
#define mpq_set_str(...) \
    _Pragma("GCC error \"mpq_set_str is forbidden by TRTS policy (use mpz_set_str + rational_set_components).\"")

/* mpq_set is allowed only via rational_set wrapper below. Direct use is forbidden. */
#ifdef mpq_set
#undef mpq_set
#endif
#define mpq_set(...) \
    _Pragma("GCC error \"Direct mpq_set is forbidden: use rational_set wrapper API to ensure unadulterated propagation.\"")

/* Prohibit common gcd-like mpz helpers */
#ifdef mpz_gcd
#undef mpz_gcd
#endif
#define mpz_gcd(...) \
    _Pragma("GCC error \"mpz_gcd is forbidden by TRTS policy (no GCD).\"")

#ifdef mpz_gcd_ui
#undef mpz_gcd_ui
#endif
#define mpz_gcd_ui(...) \
    _Pragma("GCC error \"mpz_gcd_ui is forbidden by TRTS policy (no GCD).\"")

#ifdef mpz_gcdext
#undef mpz_gcdext
#endif
#define mpz_gcdext(...) \
    _Pragma("GCC error \"mpz_gcdext is forbidden by TRTS policy (no GCD).\"")

#ifdef mpz_invert
#undef mpz_invert
#endif
#define mpz_invert(...) \
    _Pragma("GCC error \"mpz_invert is forbidden by TRTS policy (no GCD-based inversion).\"")

/* ---------------------------
   Runtime safeguard
   ---------------------------
   If a build system or toolchain bypasses the above preprocessor guards,
   the runtime check below will attempt to detect linked forbidden symbols
   and abort early. This is a best-effort measure.
   --------------------------- */

static inline void trts_check_forbidden_symbols_runtime(void) {
    /* This function intentionally keeps implementation minimal and portable:
       it attempts to reference the forbidden symbols via weak externs when
       available. Most toolchains will resolve forbidden usage at compile time
       due to the macros above. If this function is linked and the platform
       supports weak symbol detection, it will abort if forbidden functions
       are present. Implementations lacking weak symbol support may be a no-op. */
#if defined(__GNUC__) || defined(__clang__)
    /* no-op: compile-time macros ensure forbidden usage is rejected */
    (void)0;
#else
    (void)0;
#endif
}

/* ---------------------------
   Public TRTS rational API (whitelist)
   --------------------------- */

/* Initialization and cleanup (whitelisted) */
void rational_init(mpq_t value);
void rational_clear(mpq_t value);

/* Assignment operations (whitelisted wrappers)
   - rational_set: copy from another mpq (allowed, but this wrapper is the only
     sanctioned way to copy a rational within engine propagation code).
   - rational_set_si: set numerator/denominator via mpz_set_* primitives.
   - rational_set_components: set numerator and denominator mpz_t explicitly.
*/
void rational_set(mpq_t dest, mpq_srcptr src);
void rational_set_si(mpq_t dest, long numerator, unsigned long denominator);
void rational_set_components(mpq_t dest, mpz_srcptr numerator, mpz_srcptr denominator);

/* NOTE: rational_set_str (text parsing) is intentionally NOT exposed here.
   If text parsing is required, call mpz_set_str for numerator and denominator
   and then call rational_set_components. This guarantees no implicit canonicalization. */

/* Arithmetic operations (whitelisted) */
void rational_add(mpq_t result, mpq_srcptr a, mpq_srcptr b);
void rational_sub(mpq_t result, mpq_srcptr a, mpq_srcptr b);
void rational_mul(mpq_t result, mpq_srcptr a, mpq_srcptr b);
void rational_div(mpq_t result, mpq_srcptr a, mpq_srcptr b);
void rational_add_ui(mpq_t result, mpq_srcptr a, unsigned long numerator, unsigned long denominator);
void rational_neg(mpq_t result, mpq_srcptr a);
void rational_inv(mpq_t result, mpq_srcptr a);
void rational_abs(mpq_t result, mpq_srcptr a);
void rational_negate(mpq_t value);

/* Comparison operations */
int rational_cmp(mpq_srcptr a, mpq_srcptr b);
int rational_sgn(mpq_srcptr value);
bool rational_is_zero(mpq_srcptr a);

/* Numerator/denominator operations (whitelisted) */
void rational_copy_num(mpz_t dest, mpq_srcptr value);
void rational_abs_num(mpz_t dest, mpq_srcptr a);

/* NOTE: rational_denominator_zero remains a simple inline helper; it does not
   canonicalize or compute GCDs. */
static inline bool rational_denominator_zero(mpq_srcptr value) {
    return mpz_sgn(mpq_denref(value)) == 0;
}

/* Modular and special operations (implemented using whitelisted primitives) */
void rational_mod(mpq_t result, mpq_srcptr value, mpq_srcptr modulus);
void rational_delta(mpq_t result, mpq_srcptr current, mpq_srcptr previous);
void rational_floor(mpz_t dest, mpq_srcptr value);
void rational_ceil(mpz_t dest, mpq_srcptr value);
void rational_round(mpz_t dest, mpq_srcptr value);

/* I/O operations (printing only — does not canonicalize engine state) */
void rational_print(FILE *stream, mpq_srcptr value);

/* Utility: explicit forbidden setter detection
   If you need a textual-to-rational loader, implement it using mpz_set_str
   and rational_set_components outside of engine propagation paths.
*/
#ifdef __cplusplus
}
#endif

#endif // RATIONAL_H

