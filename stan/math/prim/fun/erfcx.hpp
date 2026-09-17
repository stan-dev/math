#ifndef STAN_MATH_PRIM_FUN_ERFCX_HPP
#define STAN_MATH_PRIM_FUN_ERFCX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

/**
 * Cody (1969) third-interval rational, valid for `x >= 4`.
 *
 * Gives `erfcx` directly. `x * x` is infinite for `x` large enough, which
 * correctly collapses this to the leading term `INV_SQRT_PI / x` and, at
 * infinity, to zero.
 *
 * @param x argument, `x >= 4`
 * @return scaled complementary error function
 */
inline double erfcx_cody_tail(double x) {
  const double u = 1.0 / (x * x);
  double p = 0.0163153871373020978498;
  p = 0.305326634961232344035 + u * p;
  p = 0.360344899949804439429 + u * p;
  p = 0.125781726111229246204 + u * p;
  p = 0.0160837851487422766278 + u * p;
  p = 0.000658749161529837803157 + u * p;
  double q = -1.0;
  q = -2.56852019228982242072 + u * q;
  q = -1.87295284992346047209 + u * q;
  q = -0.527905102951428412248 + u * q;
  q = -0.0605183413124413191178 + u * q;
  q = -0.00233520497626869185443 + u * q;
  return (INV_SQRT_PI + (p / q) * u) / x;
}

/**
 * Cody (1969) second-interval rational, valid for `0.46875 <= x <= 4`.
 *
 * Cody writes this interval as `erfc(x) = exp(-x*x) * R(x)`, so `erfcx`
 * is `R(x)` with the exponential cancelled analytically: no `exp` and no
 * `erfc` call. Measured at 3.7 to 7.0 ulp, of which only 0.28 ulp is
 * approximation error; the rest is rounding in the two degree-8 Horner
 * chains and the final division.
 *
 * @param y argument, `0.46875 <= y <= 4`
 * @return scaled complementary error function
 */
inline double erfcx_cody_middle(double y) {
  double p = 2.15311535474403846e-8 * y;
  p = (p + 5.64188496988670089e-1) * y;
  p = (p + 8.88314979438837594) * y;
  p = (p + 66.1191906371416295) * y;
  p = (p + 298.635138197400131) * y;
  p = (p + 881.952221241769090) * y;
  p = (p + 1712.04761263407058) * y;
  p = (p + 2051.07837782607147) * y;
  double q = y;
  q = (q + 15.7449261107098347) * y;
  q = (q + 117.693950891312499) * y;
  q = (q + 537.181101862009858) * y;
  q = (q + 1621.38957456669019) * y;
  q = (q + 3290.79923573345963) * y;
  q = (q + 4362.61909014324716) * y;
  q = (q + 3439.36767414372164) * y;
  return (p + 1230.33935479799725) / (q + 1230.33935480374942);
}

}  // namespace internal

/**
 * Return the scaled complementary error function of the argument.
 *
 * \f$\mbox{erfcx}(x) = \exp(x^2)\,\mbox{erfc}(x)\f$
 *
 * The scaling removes the Gaussian factor, so `erfcx` stays in a narrow
 * range (it decays like \f$1/(x\sqrt{\pi})\f$ for large positive `x`) where
 * the two factors individually overflow and underflow: `erfc(x)` underflows
 * to zero above `x = 27` while `exp(x*x)` overflows above `x = 26.6`, so
 * their product cannot be formed directly in the upper tail at all. It is
 * the natural building block for normal tail quantities, e.g. the standard
 * normal log CDF is
 * `LOG_HALF + log(erfcx(-x * INV_SQRT_TWO)) - x * x / 2` and the inverse
 * Mills ratio is `SQRT_TWO_OVER_SQRT_PI / erfcx(-x * INV_SQRT_TWO)`, both
 * without cancellation.
 *
 * Three branches. `erfc` dominates the cost of forming
 * `exp(x*x) * erfc(x)`, so each interval that can be served by a direct
 * approximation of `erfcx` instead is much faster.
 *
 * - `x >= 4`: Cody (1969) third-interval rational.
 * - `0.46875 <= x < 4`: Cody (1969) second-interval rational. Cody writes
 *   this interval as `erfc(x) = exp(-x*x) * R(x)`, so `erfcx` is `R(x)`
 *   and the exponential is cancelled analytically rather than computed and
 *   divided out. No `exp` and no `erfc` call.
 * - `x < 0.46875`: `exp(x * x) * erfc(x)`, with the rounding error of
 *   `x * x` recovered by `fma` and folded back in, because `exp`
 *   amplifies it into roughly `x * x * eps` -- 512 ulp at `x = -26` if
 *   left alone.
 *
 * W. J. Cody, Math. Comp. 23(107):631-638 (1969).
 * https://doi.org/10.1090/S0025-5718-1969-0247736-4
 *
 * The crossover at 4 is the same one R's `pnorm` uses: it switches to the
 * Cody tail form at `|y| > M_SQRT_32`, and since the argument here is the
 * `erfc` argument `y / sqrt(2)`, that is `sqrt(32) / sqrt(2) = 4` exactly.
 * Evaluated at 60 significant digits, this coefficient set holds to 8.6e-17
 * relative down to 4 and degrades past about 3.5, so 4 sits just inside its
 * range.
 *
 * Measured worst case over the whole range is 7.0 ulp.
 *
 * The derivative is
 * \f$\frac{d}{dx}\mbox{erfcx}(x) = 2x\,\mbox{erfcx}(x) -
 \frac{2}{\sqrt{\pi}}\f$.
 *
   \f[
   \mbox{erfcx}(x) =
   \begin{cases}
     \exp(x^2)\operatorname{erfc}(x) & \mbox{if } -\infty\leq x \leq \infty
     \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T An arithmetic type.
 * @param xx argument
 * @return scaled complementary error function of the argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double erfcx(T&& xx) {
  const double x = static_cast<double>(xx);
  constexpr double cody_min = 4.0;
  constexpr double middle_min = 0.46875;
  if (x >= cody_min) {
    return internal::erfcx_cody_tail(x);
  }
  if (x >= middle_min) {
    return internal::erfcx_cody_middle(x);
  }
  // erfcx(x) = 2 * exp(x * x) * (1 + o(1)) as x -> -infinity, which leaves
  // the binary64 range below -26.63. Returning here also keeps the branch
  // below from forming inf * 0 on -infinity.
  constexpr double overflow_max = -27.0;
  if (x < overflow_max) {
    return INFTY;
  }
  // fma(x, x, -h) is the exact rounding error of h = x * x, so the factor
  // (1 + that) restores what exp(h) would otherwise lose. The OpenCL
  // device function uses the same form; it cannot use a Dekker split
  // because the compiler simplifies t - (t - x) to x.
  const double h = x * x;
  return std::exp(h) * (1.0 + std::fma(x, x, -h)) * std::erfc(x);
}

/**
 * Structure to wrap `erfcx()` so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Scaled complementary error function applied to x.
 */
struct erfcx_fun {
  template <typename T>
  static inline auto fun(T&& x) {
    return erfcx(std::forward<T>(x));
  }
};

/**
 * Returns the elementwise `erfcx()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return Scaled complementary error function applied to each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_container_t<T>* = nullptr, require_not_var_matrix_t<T>* = nullptr>
inline auto erfcx(T&& x) {
  return apply_scalar_unary<erfcx_fun, T>::apply(std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
