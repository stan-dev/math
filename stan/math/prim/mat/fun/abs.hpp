#ifndef STAN_MATH_PRIM_MAT_FUN_ABS_HPP
#define STAN_MATH_PRIM_MAT_FUN_ABS_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
  namespace math {

    struct abs_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::fabs;
        return abs(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<abs_fun, T>::return_t
    abs(const T& x) {
      return apply_scalar_unary<abs_fun, T>::apply(x);
    }

  }
}

#endif
