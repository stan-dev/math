#ifndef STAN_MATH_PRIM_SCAL_FUN_AS_BOOL_HPP
#define STAN_MATH_PRIM_SCAL_FUN_AS_BOOL_HPP

namespace stan {
  namespace math {

    /**
     * Return true if the argument is unequal to zero and false otherwise.
     *
     * @tparam type of scalar
     * @param x value
     * @return true if value is equal to zero or NaN
     */
    template <typename T>
    inline bool as_bool(const T& x) {
      return x != 0;
    }

  }
}

#endif
