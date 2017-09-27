#ifndef STAN_MATH_PRIM_SCAL_FUN_SIZE_ZERO_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SIZE_ZERO_HPP

namespace stan {
  namespace math {

    /**
     * Returns 1 if any inputs are of length 0, returns 0
     * otherwise
     *
     * @param x arguments
     * @return 0 or 1
     */
    template <typename T>
    inline bool size_zero(const T& x) {
      return !stan::length(x);
    }

    template <typename T, typename ... Ts>
    inline bool size_zero(T& x, Ts&& ... xs) {
      return (size_zero(x) || size_zero(std::forward<Ts>(xs)...));
    }
  }
}

#endif
