#ifndef STAN_MATH_PRIM_MAT_FUN_INV_PHI_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_PHI_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/inv_Phi.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of inv_Phi().
     */
    struct inv_Phi_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using stan::math::inv_Phi;
        stan::math::check_bounded<T, T, T>(
          "inv_Phi vectorize", "Probability variable", x, 0, 1);
        return inv_Phi(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<inv_Phi_fun, T>::return_t
    inv_Phi(const T& x) {
      return apply_scalar_unary<inv_Phi_fun, T>::apply(x);
    }

  }
}

#endif
