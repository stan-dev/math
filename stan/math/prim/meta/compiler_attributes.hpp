#ifndef STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP
#define STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP

#ifdef __GNUC__
#ifndef likely
#define likely(x) __builtin_expect(!!(x), 1)
#endif
#ifndef unlikely
#define unlikely(x) __builtin_expect(!!(x), 0)
#endif
#ifdef __has_attribute
#if __has_attribute(noinline) && __has_attribute(cold)
#ifndef STAN_COLD_PATH
/**
 * Functions tagged with this attribute are not inlined and moved
 *  to a cold branch to tell the CPU to not attempt to pre-fetch
 *  the associated function.
 */
#define STAN_COLD_PATH __attribute__((noinline, cold))
#endif
#endif
#endif
#endif
#ifndef STAN_COLD_PATH
#define STAN_COLD_PATH
#endif
#ifndef likely
#define likely(x) x
#endif
#ifndef unlikely
#define unlikely(x) x
#endif

/**
 * Turns all range and size checks into no-ops
 */
#ifndef STAN_NO_RANGE_CHECKS_RETURN
/**
 * If defined, will turn off all range and size checks.
 */
#ifdef STAN_NO_RANGE_CHECKS
#define STAN_NO_RANGE_CHECKS_RETURN return
#endif
#ifndef STAN_NO_RANGE_CHECKS_RETURN
#define STAN_NO_RANGE_CHECKS_RETURN
#endif
#endif

#endif
