#ifndef STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP
#define STAN_MATH_PRIM_META_COMPILER_ATTRIBUTES_HPP

#ifdef __GNUC__
  #define likely(x) __builtin_expect(!!(x), 1)
  #define unlikely(x) __builtin_expect(!!(x), 0)
  #ifdef __has_attribute
    #if __has_attribute(noinline) && __has_attribute(cold)
      /**
       * Functions tagged with this attribute are not inlined and moved
       *  to a cold branch to tell the CPU to not attempt to pre-fetch
       *  the associated function.
       */
      #define STAN_COLD_PATH __attribute__((noinline, cold))
    #else
      #define STAN_COLD_PATH
    #endif
  #else
    #define STAN_COLD_PATH
  #endif
#else
#define likely(x) (x)
#define unlikely(x) (x)
#define STAN_COLD_PATH
#endif

/**
 * Turns all range and size checks into no-ops
 */
#ifndef STAN_NO_RANGE_AND_SIZE_CHECK
/**
 * If defined, will turn off all range and size checks.
 */
#ifdef STAN_NDEBUG
#define STAN_NO_RANGE_AND_SIZE_CHECK return
#else
#define STAN_NO_RANGE_AND_SIZE_CHECK
#endif
#else
#define STAN_NO_RANGE_AND_SIZE_CHECK
#endif
#endif
