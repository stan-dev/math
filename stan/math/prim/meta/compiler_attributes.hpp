#ifndef STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP
#define STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP

#ifdef __GNUC__
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#define STAN_COLD_PATH __attribute__((noinline, cold))
#else
#define likely(x) (x)
#define unlikely(x) (x)
#define STAN_COLD_PATH
#endif

#endif
