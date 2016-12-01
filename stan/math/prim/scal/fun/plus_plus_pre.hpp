#ifndef STAN_MATH_PRIM_SCAL_FUN_PLUS_PLUS_PRE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_PLUS_PLUS_PRE_HPP

namespace stan {
  namespace math {

    /**
     * Increment the specified argument by one and return the
     * resulting value.  Has the same effect as <code>++x</code>.
     *
     * @tparam T argument type
     * @param x argument reference
     * @return argument after having been incremented
     */
    template <typename T>
    inline T& plus_plus_pre(T& x) {
      return ++x;
    }

  }
}
#endif
