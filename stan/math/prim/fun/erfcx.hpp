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
 * Correction factor of the Cody (1969) third-interval rational:
 * `erfcx(x) = (INV_SQRT_PI + u * correction(u)) / x` with `u = 1 / x^2`.
 *
 * Split out from the value so that the derivative can reuse it. See
 * `erfcx_derivative`.
 *
 * @tparam T scalar type
 * @param u inverse square of the argument, `0 <= u <= 1/16`
 * @return `P(u) / Q(u)`
 */
template <typename T>
inline T erfcx_tail_correction(const T& u) {
  static constexpr std::array p{
      0.000658749161529837803157, 0.0160837851487422766278,
      0.125781726111229246204, 0.360344899949804439429,
      0.305326634961232344035};
  static constexpr std::array q{
      -0.00233520497626869185443, -0.0605183413124413191178,
      -0.527905102951428412248, -1.87295284992346047209,
      -2.56852019228982242072};
  T numerator = 0.0163153871373020978498;
  T denominator = -1.0;
  for (int i = 4; i >= 0; --i) {
    numerator = numerator * u + p[i];
    denominator = denominator * u + q[i];
  }
  return numerator / denominator;
}

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
  return (INV_SQRT_PI + u * erfcx_tail_correction(u)) / x;
}

/**
 * Derivative of `erfcx`, `2 * x * erfcx(x) - 2 / sqrt(pi)`.
 *
 * That difference cancels for large `x`: both terms approach
 * `2 / sqrt(pi)` while the result decays like `1 / (sqrt(pi) * x^2)`.
 * Measured against a 50-digit reference, the difference form gives 6.1 ulp
 * at `x = 4` and 2.55e+11 ulp at `x = 1e6`.
 *
 * For `x >= 4` the tail rational gives the derivative with no subtraction,
 * because the constant cancels analytically:
 *
 *   `2 * x * (INV_SQRT_PI + u * C(u)) / x - 2 / sqrt(pi) = 2 * u * C(u)`
 *
 * since `2 * INV_SQRT_PI` is `2 / sqrt(pi)`. That form measures 0.1 to
 * 31 ulp over the same range.
 *
 * @tparam T scalar type
 * @param x argument
 * @param value `erfcx(x)`
 * @return derivative of `erfcx` at `x`
 */
template <typename T>
inline T erfcx_derivative(const T& x, const T& value) {
  if (x >= 4.0) {
    const T u = 1.0 / (x * x);
    return 2.0 * u * erfcx_tail_correction(u);
  }
  return 2.0 * x * value - TWO_OVER_SQRT_PI;
}

/**
 * Cody (1969) second-interval rational, valid for `0.46875 <= x <= 4`.
 *
 * Cody writes this interval as `erfc(x) = exp(-x*x) * R(x)`, so `erfcx`
 * is `R(x)` with the exponential cancelled analytically: no `exp` and no
 * `erfc` call. Measured at 3.7 to 7.0 ulp.
 *
 * @param y argument, `0.46875 <= y <= 4`
 * @return scaled complementary error function
 */
inline double erfcx_cody_middle(double y) {
  static constexpr std::array p{1230.33935479799725, 2051.07837782607147,
                                1712.04761263407058, 881.952221241769090,
                                298.635138197400131, 66.1191906371416295,
                                8.88314979438837594, 5.64188496988670089e-1};
  static constexpr std::array q{1230.33935480374942, 3439.36767414372164,
                                4362.61909014324716, 3290.79923573345963,
                                1621.38957456669019, 537.181101862009858,
                                117.693950891312499, 15.7449261107098347};
  // Pair the coefficients first, so the four products are independent, then
  // run a Horner chain half as long in y2. Written as arithmetic, not
  // std::fma: with the Stan Math flags FP_FAST_FMA is not defined, so an
  // explicit std::fma is a libm call. Measured on a Xeon E5-2680 v3, the
  // std::fma form costs 31.0 ns per call against 4.46 ns for this one. The
  // compiler contracts these into hardware fma wherever -march allows it.
  const double y2 = y * y;
  std::array<double, 4> p_vals;
  std::array<double, 4> q_vals;
  for (int i = 0, j = 0; i < 4; ++i, j += 2) {
    p_vals[i] = p[j] + p[j + 1] * y;
    q_vals[i] = q[j] + q[j + 1] * y;
  }
  double num = p_vals[3] + 2.15311535474403846e-8 * y2;
  double den = y2 + q_vals[3];
  for (int i = 2; i >= 0; --i) {
    num = p_vals[i] + num * y2;
    den = q_vals[i] + den * y2;
  }
  return num / den;
}

/**
 * Degree-18 polynomial for `|x| < 0.46875`.
 *
 * A Chebyshev-economized expansion of `erfcx` about zero, so it needs no
 * library call and no branch on the sign of `x`. Its low-order coefficients
 * reproduce the Maclaurin series of `erfcx` exactly (`1`, `-2/sqrt(pi)`,
 * `1`, ...), which is a useful check that the fit is right. `even[0]`,
 * `even[1]` and `odd[0]` are those three terms.
 *
 * Measured at 2.2 ulp over the interval. Evaluating the same coefficients
 * as one degree-18 Horner chain gives 1.5 ulp instead, but runs 1.7 times
 * slower, because that chain is 18 dependent operations. The whole function
 * is 7.0 ulp, set by the middle interval, so the 0.7 ulp costs nothing.
 *
 * The plain Maclaurin series is only usable to about `|x| = 0.125`; four
 * further economized terms extend it to 0.46875 at no measurable cost,
 * which is what removes the last `exp` call from the positive axis.
 *
 * @param x argument, `|x| < 0.46875`
 * @return scaled complementary error function
 */
inline double erfcx_small(double x) {
  // Split into the even and odd powers of x, so the two Horner chains run
  // independently. A single degree-18 chain is 18 dependent operations; two
  // chains of 9 halve that latency.
  static constexpr std::array even = {1.0,
                                      4.99999999999992839e-01,
                                      1.66666666667239644e-01,
                                      4.16666666458337179e-02,
                                      8.33333374332981443e-03,
                                      1.38888415444527033e-03,
                                      1.98445679338826757e-04,
                                      2.46655529768908249e-05};
  static constexpr std::array odd
      = {-1.12837916709551256,     -7.52252778063651983e-01,
         -3.00901111227312890e-01, -8.59717459974174147e-02,
         -1.91048337772546720e-02, -3.47359067853470795e-03,
         -5.34506929034156810e-04, -7.08163358203131886e-05};
  static_assert(even.size() == odd.size(),
                "the two chains are stepped by one loop");
  // The leading coefficient of each chain seeds the accumulator, and the
  // trailing 1.0 is added at the end, so both arrays hold the interior
  // coefficients only. Derive the bound from the array rather than writing
  // it out, so shortening an array cannot leave the loop reading past it.
  const double x2 = x * x;
  double e = 3.05977060678449757e-06;
  double o = -9.35890030086883823e-06;
  for (int i = static_cast<int>(even.size()) - 1; i >= 0; --i) {
    e = even[i] + x2 * e;
    o = odd[i] + x2 * o;
  }
  return 1.0 + x2 * e + x * o;
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
 * Four branches. The whole positive axis is covered without calling any
 * library function, which is what makes this fast: a direct approximation
 * of `erfcx` avoids `erfc`, and `erfc` dominates the cost of forming
 * `exp(x*x) * erfc(x)`.
 *
 * - `x >= 4`: Cody (1969) third-interval rational.
 * - `0.46875 <= x < 4`: Cody (1969) second-interval rational. Cody writes
 *   this interval as `erfc(x) = exp(-x*x) * R(x)`, so `erfcx` is `R(x)`
 *   and the exponential is cancelled analytically rather than computed and
 *   divided out.
 * - `|x| < 0.46875`: a degree-18 Chebyshev-economized expansion about
 *   zero, covering both signs with no branch and no library call.
 * - `x <= -0.46875`: the reflection `2*exp(x*x) - erfcx(-x)`, which is the
 *   one place a library call is unavoidable, since `erfcx` grows like
 *   `2*exp(x*x)` as `x -> -infinity`. The rounding error of `x * x` is
 *   recovered with `fma` and folded back in, because `exp` amplifies it
 *   into roughly `x * x * eps` -- 512 ulp at `x = -26` if left alone.
 *   Below `x = -6.1`, `erfcx(-x)` is under `eps/2` of `2*exp(x*x)`, so
 *   the subtraction is skipped.
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
  constexpr double overflow_max = -27.0;
  if (x >= cody_min) {
    return internal::erfcx_cody_tail(x);
  } else if (x >= middle_min) {
    return internal::erfcx_cody_middle(x);
  } else if (x > -middle_min) {
    return internal::erfcx_small(x);
  } else if (x < overflow_max) {
    return INFTY;
  }
  const double h = x * x;
  const double two_exp_x2 = 2.0 * std::exp(h) * (1.0 + std::fma(x, x, -h));
  constexpr double reflect_min = -6.1;
  if (x < reflect_min) {
    return two_exp_x2;
  }
  const double y = -x;
  return two_exp_x2
         - (y >= cody_min ? internal::erfcx_cody_tail(y)
                          : internal::erfcx_cody_middle(y));
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
