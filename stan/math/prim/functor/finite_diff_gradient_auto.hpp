#ifndef STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP
#define STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument using finite differences.
 *
 * @tparam F Type of function
 * @param[in] f function
 * @param[in] x argument to function
 * @param[out] fx function applied to argument
 * @param[out] grad_fx gradient of function at argument
 */
template <typename F, typename VectorT,
          typename ScalarT = return_type_t<VectorT>>
void finite_diff_gradient_auto(const F& f, const VectorT& x, ScalarT& fx,
                               VectorT& grad_fx) {
  using boost::math::differentiation::finite_difference_derivative;
  VectorT x_temp(x);
  fx = f(x);
  grad_fx.resize(x.size());

  for (Eigen::Index i = 0; i < x.size(); ++i) {
    auto fun = [&i, &x_temp, &f](const auto& y) {
      x_temp[i] = y;
      return f(x_temp);
    };
    grad_fx[i] = finite_difference_derivative(fun, x[i]);
    x_temp[i] = x[i];
  }
}

}  // namespace math
}  // namespace stan
#endif
