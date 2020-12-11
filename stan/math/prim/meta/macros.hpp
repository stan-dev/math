#ifndef STAN_MATH_PRIM_META_MACROS_HPP
#define STAN_MATH_PRIM_META_MACROS_HPP

#ifdef __GNUC__
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#define STAN_ALWAYS_INLINE __attribute__((always_inline))
#define STAN_STRONG_INLINE inline STAN_ALWAYS_INLINE
#define STAN_FLATTEN __attribute__((flatten))
#define STAN_STRONGER_INLINE inline __attribute__((always_inline, flatten))
#define STAN_COLD_PATH __attribute__((noinline, cold))
#define STAN_HOT_PATH __attribute__((hot))
#else
#define likely(x) (x)
#define unlikely(x) (x)
#define STAN_ALWAYS_INLINE
#define STAN_STRONG_INLINE inline
#define STAN_FLATTEN
#define STAN_STRONGER_INLINE
#define STAN_COLD_PATH
#define STAN_HOT_PATH
#endif


#ifdef STAN_REMOVE_RANGE_AND_SIZE_CHECKS
#define STAN_NO_RANGE_AND_SIZE_CHECK return;
#else
#define STAN_NO_RANGE_AND_SIZE_CHECK
#endif

#endif
