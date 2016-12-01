#ifndef STAN_MATH_PRIM_MAT_FUN_PLUS_PLUS_PRE_HPP
#define STAN_MATH_PRIM_MAT_FUN_PLUS_PLUS_PRE_HPP

#include <stan/math/prim/mat/fun/add.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
  namespace math {

    /**
     * Increment the elements of the specified argument by one and
     * return the resulting value.  Delegates to the <code>add</code>
     * function with second argument <code>1</code>.
     *
     * @tparam T scalar type for argument
     * @tparam R row specification for argument
     * @tparam C column specification for argument
     * @param x argument reference
     * @return argument after having been incremented
     */
    template <typename T, int R, int C>
    inline Eigen::Matrix<T, R, C>& plus_plus_pre(Eigen::Matrix<T, R, C>& x) {
      x = add(x, 1);
      return x;
    }

  }
}
#endif
