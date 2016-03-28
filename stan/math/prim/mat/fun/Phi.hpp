#ifndef STAN_MATH_PRIM_MAT_FUN_PHI_HPP
#define STAN_MATH_PRIM_MAT_FUN_PHI_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/Phi.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of Phi().
     */
    struct Phi_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::Phi;
        return Phi(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<Phi_fun, T>::return_t
    Phi(const T& x) {
      return apply_scalar_unary<Phi_fun, T>::apply(x);
    }

  }
}

#endif
