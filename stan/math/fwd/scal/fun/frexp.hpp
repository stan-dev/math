#ifndef STAN_MATH_FWD_SCAL_FUN_FREXP_HPP
#define STAN_MATH_FWD_SCAL_FUN_FREXP_HPP

#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Decomposes a float-like variable into a normalized 
     * fraction and an integral power of two. (cmath)
     *
     * @tparam T type of scalar of the variable that will be
     * decomposed.
     * @param[in] a The variable that will be decomposed.
     * @param[out] b Pointer to an integer that will stored the
     * integral power of two.
     * @return Normalized fraction
     */
      template <typename T>
      inline fvar<T> frexp(const fvar<T>& x, int* b) {
          using std::frexp;
          return fvar<T>(frexp(x.val_, b), 0);
      }

  }
}
#endif
