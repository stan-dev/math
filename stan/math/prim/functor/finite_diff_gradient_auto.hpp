#ifndef STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP
#define STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/finite_diff_stepsize.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument using finite difference.
 *
 * <p>The functor must implement
 *
 * <code>
 * double operator()(const Eigen::Matrix<double, -, 1>&) const;
 * </code>
 *
 * <p>Error of derivative in dimension `i` should be on the should be on
 * order of `epsilon(i)^6`, where `epsilon(i) = sqrt(delta) * abs(x(i))`
 * for input `x` at dimension `i`.
 *
 * The reference for this algorithm is:
 *
 * <br />Robert de Levie. 2009. An improved numerical approximation
 * for the first derivative. Journal of Chemical Sciences 121(5), page
 * 3.
 *
 * <p>The reference for automatically setting the difference is this
 * section of the Wikipedia,
 *
 * <br /><a
 * href="https://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating-point_arithmetic">Numerical
 * differentiation: practical considerations using floating point
 * arithmetic</a>.
 *
 * <p>Evaluating this function involves 6 calls to the function being
 * differentiated for each dimension in the input, plus one global
 * evaluation.  All evaluations will be for double-precision inputs.
 *
 * @tparam F Type of function
 * @param[in] f function
 * @param[in] x argument to function
 * @param[out] fx function applied to argument
 * @param[out] grad_fx gradient of function at argument
 */
template <typename F, typename VectorT, typename GradVectorT,
          typename ScalarT = return_type_t<VectorT>>
void finite_diff_gradient_auto(const F& f, VectorT&& x, ScalarT& fx,
                               GradVectorT& grad_fx) {
  using EigT = Eigen::Matrix<ScalarT, -1, 1>;
  static constexpr int h_scale[6] = {3, 2, 1, -3, -2, -1};
  static constexpr int mults[6] = {1, -9, 45, -1, 9, -45};

  fx = f(x);
  grad_fx.resize(x.size());
  Eigen::Map<EigT> grad_map(grad_fx.data(), grad_fx.size());

  grad_map = EigT::NullaryExpr(x.size(), [&f, &x](Eigen::Index i) {
    double h = finite_diff_stepsize(value_of_rec(x[i]));
    ScalarT delta_f = 0;
    for (int j = 0; j < 6; ++j) {
      auto x_temp
          = EigT::NullaryExpr(x.size(), [&x, &i, &h, &j](Eigen::Index k) {
              return k == i ? x[i] + h * h_scale[j] : x[k];
            });
      delta_f += f(std::move(x_temp)) * mults[j];
    }
    return delta_f / (60 * h);
  });
}

}  // namespace math
}  // namespace stan
#endif
