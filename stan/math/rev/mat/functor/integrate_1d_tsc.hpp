#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_integrate_1d_tsc_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_integrate_1d_tsc_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/to_var.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/rev/arr/fun/to_var.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>

#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/rev/mat/functor/de_integrator.hpp>
#include <functional>
#include <cmath>
#include <ostream>
#include <vector>

namespace stan {

namespace math {

/**
 * Return the numerical integral of a function f given its gradient
 * g.
 *
 * This function uses the algorithm of integrate_1d_tscg for
 * numerical integration.
 *
 * @tparam T Type of f.
 * @tparam G Type of g.
 * @param f a functor with signature
 * double (double, std::vector<T_param>, std::vector<T_x_r>,
 * std::vector<T_x_i>, std::ostream&) or with signature
 * double (double, T_param, T_x_r, T_x_i, std::ostream&)
 * or a "mixture" of these, where the first argument is one being
 * integrated and the second one is either an extra scalar or vector
 * being passed to f.
 * @param g a functor with signature
 * double (double, std::vector<T_param>, std::vector<T_x_r>,
 * std::vector<T_x_i>, int, std::ostream&) or with signature double
 * (double, T_param, T_x_r, T_x_i, int, std::ostream&)
 * where the first argument is being integrated and the second one is
 * either an extra scalar or vector being passed to f and the
 * fifth one selects which component of the gradient vector
 * is to be returned.
 * @param a lower limit of integration, must be a finite double.
 * @param b upper limit of integration, must be a finite double.
 * @param tre target relative error.
 * @param tae target absolute error.
 * @param param additional parameters to be passed to f.
 * @param T_x_r additional data to be passed to f.
 * @param T_x_i additional integer data to be passed to f.
 * @param msgs stream.
 * @return numeric integral of function f.
 */
template <typename F, typename G, typename T_param, typename T_x_r,
          typename T_x_i>
inline typename scalar_type<T_param>::type integrate_1d_tscg(
    const F& f, const G& g, const double a, const double b,
    const T_param& param, const T_x_r& x_r, const T_x_i& x_i,
    std::ostream& msgs, const double tre = 1e-6, const double tae = 1e-6) {
  check_finite("integrate_1d_tsc", "lower limit", a);
  check_finite("integrate_1d_tsc", "upper limit", b);

  using std::placeholders::_1;

  // hard case, we want a normalizing factor
  if (!is_constant_struct<T_param>::value) {
    size_t N = length(param);
    std::vector<double> grad(N);

    auto value_of_param = value_of(param);

    for (size_t n = 0; n < N; n++)
      grad[n] = de_integrator(
          std::bind<double>(g, _1, value_of_param, x_r, x_i,
                            static_cast<int>(n + 1), std::ref(msgs)),
          a, b, tre, tae);

    double val_ = de_integrator(
        std::bind<double>(f, _1, value_of_param, x_r, x_i, std::ref(msgs)), a,
        b, tre, tae);

    operands_and_partials<T_param> ops_partials(param);
    for (size_t n = 0; n < N; n++)
      ops_partials.edge1_.partials_[n] += grad[n];

    return ops_partials.build(val_);
  }
  // easy case, here we are calculating a normalizing constant,
  // not a normalizing factor, so g doesn't matter at all
  return de_integrator(
      std::bind<double>(f, _1, value_of(param), x_r, x_i, std::ref(msgs)), a, b,
      tre, tae);
}

/**
 * Calculate gradient of f(x, param, std::ostream&)
 * with respect to param_n (which must be an element of param)
 */
template <typename F, typename T_param, typename T_x_r, typename T_x_i>
inline double gradient_of_f(const F& f, const double x, const T_param& param,
                            const T_x_r& x_r, const T_x_i& x_i,
                            const var& param_n, std::ostream& msgs) {
  set_zero_all_adjoints_nested();
  f(x, param, x_r, x_i, msgs).grad();
  return param_n.adj();
}

/**
 * Return the numerical integral of a function f with its
 * gradients being inferred automatically (but slowly).
 *
 * The numerical integration algorithm used is the Tanh-sinh
 * quadrature method (also known as Double Exponential
 * Transformation) which was proposed by Hidetosi Takahasi and
 * Masatake Mori in 1974
 *
 * The implementation of integration used is given John D. Cook with
 * the slight modification of adding a relative tolerance error. See
 * www.codeproject.com/kb/recipes/fastnumericalintegration.aspx
 * www.johndcook.com/blog/double_exponential_integration/
 *
 * Such implementation assumes f to be smooth with no discontinuity
 * in the function nor in any of its derivatives.
 *
 * The integration is terminated when both relative and absolute
 * tolerance are reached or when a predefined number of iterations
 * is reached (such termination condition was designed to target
 * floating point).
 *
 * @tparam T Type of f.
 * @param f a functor with signature
 * double (double, std::vector<T_param>, std::vector<T_x_r>,
 * std::vector<T_x_i>, std::ostream&) or with signature
 * double (double, T_param, T_x_r, T_x_i, std::ostream&)
 * or a "mixture" of these, where the first argument is one being
 * integrated and the second one is either an extra scalar or vector
 * being passed to f.
 * @param a lower limit of integration, must be double type.
 * @param b upper limit of integration, must be double type.
 * @param tre target relative error.
 * @param tae target absolute error.
 * @param param additional parameters to be passed to f.
 * @param T_x_r additional data to be passed to f.
 * @param T_x_i additional integer data to be passed to f.
 * @param msgs stream.
 * @return numeric integral of function f.
 */
template <typename F, typename T_param, typename T_x_r, typename T_x_i>
inline typename scalar_type<T_param>::type integrate_1d_tsc(
    const F& f, const double a, const double b, const T_param& param,
    const T_x_r& x_r, const T_x_i& x_i, std::ostream& msgs,
    const double tre = 1e-6, const double tae = 1e-6) {
  using std::placeholders::_1;

  stan::math::check_finite("integrate_1d_tsc", "lower limit", a);
  stan::math::check_finite("integrate_1d_tsc", "upper limit", b);

  double val_ = de_integrator(
      std::bind<double>(f, _1, value_of(param), x_r, x_i, std::ref(msgs)), a, b,
      tre, tae);

  if (!is_constant_struct<T_param>::value) {
    size_t N = stan::length(param);
    std::vector<double> results(N);

    try {
      start_nested();

      auto clean_param = to_var(value_of(param));
      typedef decltype(clean_param) clean_T_param;

      scalar_seq_view<clean_T_param> clean_param_vec(clean_param);

      for (size_t n = 0; n < N; n++)
        results[n] = de_integrator(
            std::bind<double>(gradient_of_f<F, clean_T_param, T_x_r, T_x_i>, f,
                              _1, clean_param, x_r, x_i, clean_param_vec[n],
                              std::ref(msgs)),
            a, b, tre, tae);
    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();

    operands_and_partials<T_param> ops_partials(param);
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
