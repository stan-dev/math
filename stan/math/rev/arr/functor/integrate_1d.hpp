#ifndef STAN_MATH_PRIM_REV_FUNCTOR_integrate_1d_HPP
#define STAN_MATH_PRIM_REV_FUNCTOR_integrate_1d_HPP

#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/rev/arr/fun/to_var.hpp>
#include <stan/math/prim/scal/functor/de_integrator.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/to_var.hpp>
#include <functional>
#include <cmath>
#include <ostream>
#include <vector>

namespace stan {

namespace math {

/**
 * Calculate gradient of f(x, param, std::ostream&)
 * with respect to param_n (which must be an element of param)
 */
template <typename F>
inline double gradient_of_f(const F& f, const double x,
                            const std::vector<var>& param,
                            const std::vector<double>& x_r,
                            const std::vector<int>& x_i, const var& param_n,
                            std::ostream& msgs) {
  set_zero_all_adjoints_nested();
  f(x, param, x_r, x_i, msgs).grad();
  return param_n.adj();
}

/**
 * Return the numerical integral of a function f.
 *
 * The limits of integration cannot be autodiff variables.
 *
 * The numerical integration algorithm used is the Tanh-sinh
 * quadrature method (also known as Double Exponential
 * Transformation) which was proposed by Hidetosi Takahasi and
 * Masatake Mori in 1974.
 *
 * The implementation of integration used is given John D. Cook. See
 * www.codeproject.com/kb/recipes/fastnumericalintegration.aspx
 * www.johndcook.com/blog/double_exponential_integration/
 *
 * The signature for the function to be integrated is:
 * double (double x, std::vector<var> params, std::vector<double> x_r,
 *   std::vector<int> x_i, std::ostream& msgs)
 *
 * It should return the value of the function evaluated at x. Any errors should
 * be printed to the msgs stream.
 *
 * Such implementation assumes f to be smooth with no discontinuity
 * in the function nor in any of its derivatives
 *
 * The integration is terminated when the absolute value of the current estimate
 * of the integral minus the previous estimate of the integral is less than
 * tolerance. The same tolerance is applied for the gradients as well.
 *
 * @tparam T Type of f
 * @param f a functor to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param tolerance target absolute error
 * @param param additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param msgs stream
 * @return numeric integral of function f
 */
template <typename F, typename var>
inline var integrate_1d(const F& f, const double a, const double b,
                        const std::vector<var>& param,
                        const std::vector<double>& x_r,
                        const std::vector<int>& x_i, std::ostream& msgs,
                        const double tolerance = 1e-6) {
  using std::placeholders::_1;
  static const char* function = "integrate_1d";

  stan::math::check_finite(function, "lower limit", a);
  stan::math::check_finite(function, "upper limit", b);
  stan::math::check_less_or_equal(function, "lower limit", a, b);

  double val_ = de_integrator(
      std::bind<double>(f, _1, value_of(param), x_r, x_i, std::ref(msgs)), a, b,
      tolerance);

  if (!is_constant_struct<var>::value) {
    size_t N = stan::length(param);
    std::vector<double> results(N);

    try {
      start_nested();

      auto clean_param = to_var(value_of(param));

      for (size_t n = 0; n < N; n++)
        results[n] = de_integrator(
            std::bind<double>(gradient_of_f<F>, f, _1, clean_param, x_r, x_i,
                              clean_param[n], std::ref(msgs)),
            a, b, tolerance);
    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();

    operands_and_partials<std::vector<var> > ops_partials(param);
    for (size_t n = 0; n < N; n++)
      ops_partials.edge1_.partials_[n] += results[n];

    return ops_partials.build(val_);
  } else {
    return val_;
  }
}

}  // namespace math

}  // namespace stan

#endif
