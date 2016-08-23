#ifndef STAN_MATH_FWD_SCAL_FUN_FREXP_HPP
#define STAN_MATH_FWD_SCAL_FUN_FREXP_HPP

#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
  namespace math {
      
      template <typename T>
      inline fvar<T> frexp(const fvar<T>& x, int* b) {
          using std::frexp;
          return fvar<T>(frexp(x.val_, b), 0);
      }
      
  }
}
#endif
