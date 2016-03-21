#ifndef STAN_MATH_PRIM_MAT_FUN_SIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_SIN_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    struct sin_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::sin;
        return sin(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<sin_fun, T>::return_t
    sin(const T& x) {
      return apply_scalar_unary<sin_fun, T>::apply(x);
    }

  }
}

#endif
