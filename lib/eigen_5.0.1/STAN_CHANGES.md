# Stan changes to vendored Eigen 5.0.1

This directory contains Eigen 5.0.1 plus the following Stan-local patches.

## Core fill initialization

- `Eigen/src/Core/Fill.h`
- Avoid the optimized `memset` fill path when `NumTraits<Xpr>::RequireInitialization`
  is true.
- This keeps Eigen from bypassing required scalar initialization for Stan scalar
  types that are trivially copyable but still require construction semantics.

## GCC intrinsic compatibility helpers

- `Eigen/src/Core/arch/AVX/TypeCasting.h`
- `Eigen/src/Core/arch/AVX512/PacketMath.h`
- `Eigen/src/Core/arch/AVX512/TypeCasting.h`
- Route selected AVX and AVX512 intrinsics through Eigen-local helper functions
  so older GCC versions use compatible intrinsic forms.
- The guarded helpers cover `_mm256_set_m128*`, AVX512 unaligned integer
  load/store intrinsics, and `_mm512_cmpneq_ps_mask`.
- The version guards use narrower cutoffs for the affected GCC intrinsic
  availability: GCC < 8.1 for the selected AVX casts and AVX512 compare helper,
  and GCC < 10.1 for the selected AVX512 integer load/store helpers.
