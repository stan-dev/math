#ifndef STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP
#define STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP

#ifdef __GNUC__
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#ifdef __has_attribute
#if __has_attribute (noinline) && __has_attribute (cold)
#define STAN_COLD_PATH __attribute__((noinline, cold))
#endif
#endif
#else
#define likely(x) (x)
#define unlikely(x) (x)
#define STAN_COLD_PATH
#endif


#ifndef STAN_NO_RANGE_AND_SIZE_CHECK
#ifdef STAN_NDEBUG
#define STAN_NO_RANGE_AND_SIZE_CHECK return
#endif
#else
#define STAN_NO_RANGE_AND_SIZE_CHECK
#endif

#endif
