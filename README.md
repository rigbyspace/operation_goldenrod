# operation_goldenrod
Some things to satisfy my own curiosity and nothing more

# TRTS Engine - Complete Reconstruction

## Overview

This is a complete, clean reconstruction of the TRTS (Temporal Rational Transformation System) engine based on the architectural blueprint from the audit document. The engine implements deterministic, strictly rational propagation with no canonicalization, no GCD reduction, and no floating-point approximation during computation.

## Architecture

The reconstruction follows the blueprint exactly with these core modules:

### Core Modules

1. **rational.h/c** - Strict rational arithmetic
   - Zero-based counting (0/0 representation)
   - No canonicalization or simplification
   - All operations preserve raw numerator/denominator values

2. **state.h/c** - Complete TRTS state structure
   - Primary registers: υ (upsilon), β (beta), κ (koppa)
   - Supplementary registers: ε (epsilon), φ (phi)
   - 4-level koppa stack for multi-level mode
   - All flags and counters for deterministic execution

3. **config.h/c** - Configuration system
   - All mode enumerations (engine, psi, koppa)
   - Feature flags for optional behaviors
   - Rational seeds and parameters

4. **psi.h/c** - Psi transforms
   - Standard 2-way: (υ,β) → (β/υ, υ/β)
   - Triple 3-way: (υ,β,κ) → (β/κ, κ/υ, κ/β)
   - Explicit zero-denominator checks

5. **koppa.h/c** - Koppa operations
   - DUMP, POP, ACCUMULATE modes
   - 4-level FIFO stack management
   - Sample selection at specific microticks

6. **engine.h/c** - Core propagation engine
   - Three track modes: ADD, MULTI, SLIDE
   - Asymmetric cascade
   - Stack depth and koppa magnitude modulation
   - Delta cross-propagation
   - Modular wrap

7. **simulate.h/c** - Simulation orchestrator
   - 11 microticks per tick (E, M, R phases)
   - Pattern detection (primes, Fibonacci, perfect powers)
   - Ratio triggers
   - CSV output or observer callback

### Additional Modules

8. **config_loader.h/c** - JSON configuration parser
9. **analysis_utils.h/c** - In-memory statistical analysis

## TRTS Axioms (Enforced Throughout)

1. **No Canonicalization** - Rationals are never reduced or simplified
2. **No GCD Operations** - Fractions remain in original form
3. **Zero Invariant** - Zero numerators force zero denominators (0/0)
4. **Strict Rational Arithmetic** - All computation in ℚ, no floating point
5. **Deterministic Propagation** - Same inputs always produce same outputs
6. **Evaluation Outside Propagation** - Pattern checks never alter computation

## Building

### Requirements
- GCC or compatible C compiler
- GMP library (GNU Multiple Precision Arithmetic Library)
- Make

### Compilation

```bash
make all
```

This builds:
- `libtrts.a` - Core library
- `trts_simulate` - Full simulation with CSV output
- `trts_go_time` - Minimal CLI runner

### Running

**Full simulation:**
```bash
./trts_simulate --ticks 100 --ups 3/2 --beta 5/3 --triple-psi
```

**Minimal simulation:**
```bash
./trts_go_time --ticks 50 --ups 1/1 --beta 1/1 --output run.csv
```

## Microtick Sequencing

Each tick consists of 11 microticks with specific phases:

- **E (Epsilon)**: Microticks 1, 4, 7, 10
  - Compute ε = υ
  - Run engine step
  - Check for patterns
  - Microtick 10 forces ρ (emission event)

- **M (Memory)**: Microticks 2, 5, 8, 11
  - Check β for patterns
  - Decide ψ firing
  - Accrue κ
  - Reset ρ latch

- **R (Reset)**: Microticks 3, 6, 9
  - Accrue κ without ψ
  - Clear flags

## Configuration

Configurations can be loaded from JSON files:

```json
{
  "tick_count": 100,
  "engine_mode": 0,
  "psi_mode": 0,
  "triple_psi": true,
  "multi_level_koppa": true,
  "upsilon_seed": "3/2",
  "beta_seed": "5/3",
  "koppa_seed": "1/1"
}
```

## Output Files

**events.csv** - Event flags per microtick:
- tick, microtick, phase
- rho_event, psi_fired, mu_zero
- forced_emission, ratio_triggered
- stack depth, modes

**values.csv** - Raw rational values per microtick:
- tick, microtick
- upsilon_num, upsilon_den
- beta_num, beta_den
- koppa_num, koppa_den
- stack contents
- delta values
- triangle ratios

## Analysis

The analysis_utils module provides in-memory statistical analysis:
- Convergence detection
- Pattern classification
- Ratio statistics
- Psi spacing analysis
- Stack depth distribution

All analysis uses an observer callback to avoid file I/O and maintain determinism.

## Compliance

This reconstruction strictly follows the blueprint with:
- ✓ Complete rational API with 0/0 invariant
- ✓ All state registers and flags
- ✓ Complete configuration system
- ✓ Standard and triple psi transforms
- ✓ Koppa operations with 4-level stack
- ✓ Engine with all modulation features
- ✓ 11-microtick simulation loop
- ✓ Pattern detection outside propagation
- ✓ Ratio triggers
- ✓ CSV output and observer callbacks
- ✓ In-memory analysis
- ✓ No canonicalization anywhere
- ✓ Deterministic behavior

## License

GPL v3 (as specified in LICENSE)

## References

This reconstruction is based on the complete audit and architectural blueprint of the operation_goldenrod repository, implementing all modules, data structures, algorithms, invariants, and interactions as documented.
