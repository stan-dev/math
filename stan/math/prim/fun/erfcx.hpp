#ifndef STAN_MATH_PRIM_FUN_ERFCX_HPP
#define STAN_MATH_PRIM_FUN_ERFCX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>

namespace stan {
namespace math {

namespace internal {

// Cody (1969) third interval, ascending powers of `u = 1 / x^2`.
inline constexpr std::array<double, 6> erfcx_tail_p{
    0.000658749161529837803157, 0.0160837851487422766278,
    0.125781726111229246204,    0.360344899949804439429,
    0.305326634961232344035,    0.0163153871373020978498};
inline constexpr std::array<double, 6> erfcx_tail_q{
    -0.00233520497626869185443, -0.0605183413124413191178,
    -0.527905102951428412248,   -1.87295284992346047209,
    -2.56852019228982242072,    -1.0};

// above this the correction is below eps / 8 and x^12 overflows
inline constexpr double erfcx_tail_leading_only = 0x1p27;

// Cody (1969) second interval, ascending powers of `y`.
inline constexpr std::array<double, 9> erfcx_middle_p{
    1230.33935479799725, 2051.07837782607147,    1712.04761263407058,
    881.952221241769090, 298.635138197400131,    66.1191906371416295,
    8.88314979438837594, 5.64188496988670089e-1, 2.15311535474403846e-8};
inline constexpr std::array<double, 9> erfcx_middle_q{
    1230.33935480374942, 3439.36767414372164, 4362.61909014324716,
    3290.79923573345963, 1621.38957456669019, 537.181101862009858,
    117.693950891312499, 15.7449261107098347, 1.0};

/**
 * The two tail polynomials, taking the coefficients in the order the
 * iterators give them as descending powers of `t`.
 *
 * @tparam N number of coefficients in each array
 * @tparam T scalar type
 * @tparam It coefficient iterator
 * @param p first numerator coefficient
 * @param q first denominator coefficient
 * @param t argument
 * @return numerator and denominator
 */
template <std::size_t N, typename T, typename It>
inline std::pair<T, T> erfcx_tail_polynomials(It p, It q, const T& t) {
  T numerator = *p;
  T denominator = *q;
  for (std::size_t i = 1; i < N; ++i) {
    numerator = numerator * t + *++p;
    denominator = denominator * t + *++q;
  }
  return {numerator, denominator};
}

/**
 * Horner in `y^2` over adjacent coefficient pairs, so the two halves of the
 * chain are independent.
 *
 * @tparam T scalar type
 * @tparam N number of coefficients
 * @param c coefficients, ascending
 * @param y argument
 * @param y2 `y * y`
 * @return polynomial value
 */
template <typename T, std::size_t N>
inline T erfcx_paired_horner(const std::array<double, N>& c, const T& y,
                             const T& y2) {
  int j = static_cast<int>(N) - (N % 2 ? 3 : 4);
  T r = N % 2 ? T(c[N - 1]) : T(c[N - 2] + c[N - 1] * y);
  for (; j >= 0; j -= 2) {
    r = r * y2 + (c[j] + c[j + 1] * y);
  }
  return r;
}

/**
 * Cody (1969) third-interval rational, valid for `x >= 4`.
 *
 * `erfcx(x) = (INV_SQRT_PI + u * C(u)) / x` with `u = 1 / x^2` and `C` the
 * Cody correction. `u * C(u)` is formed from the coefficients in `x^2`, so
 * `1 / x^2` is never formed.
 *
 * @param x argument, `x >= 4`
 * @return scaled complementary error function
 */
inline double erfcx_cody_tail(double x) {
  if (x > erfcx_tail_leading_only) {
    return INV_SQRT_PI / x;
  }
  const double s = x * x;
  const auto pq = erfcx_tail_polynomials<erfcx_tail_p.size()>(
      erfcx_tail_p.cbegin(), erfcx_tail_q.cbegin(), s);
  return (INV_SQRT_PI + pq.first / (s * pq.second)) / x;
}

/**
 * Derivative of `erfcx`, `2 * x * erfcx(x) - 2 / sqrt(pi)`.
 *
 * That difference cancels for large `x`: both terms approach
 * `2 / sqrt(pi)` while the result decays like `1 / (sqrt(pi) * x^2)`.
 * For `x >= 4` the tail rational gives the derivative with no subtraction,
 * because the constant cancels analytically:
 *
 *   `2 * x * (INV_SQRT_PI + u * C(u)) / x - 2 / sqrt(pi) = 2 * u * C(u)`
 *
 * since `2 * INV_SQRT_PI` is `2 / sqrt(pi)`. The walk in `x^2` is the one
 * `erfcx_cody_tail` uses.
 *
 * Above 30 the derivative is taken from the asymptotic expansion of
 * `-sqrt(pi) * x^2 * erfcx'(x)` in `u = 1 / x^2`, whose coefficients
 * satisfy `c_{n+1} = -(n + 3/2) * c_n`.
 *
 * @tparam T scalar type
 * @param x argument
 * @param value `erfcx(x)`
 * @return derivative of `erfcx` at `x`
 */
template <typename T>
inline T erfcx_derivative(const T& x, const T& value) {
  if (x < 4.0) {
    return 2.0 * x * value - TWO_OVER_SQRT_PI;
  }
  const T x2 = x * x;
  if (x < 30.0) {
    const auto pq = erfcx_tail_polynomials<erfcx_tail_p.size()>(
        erfcx_tail_p.cbegin(), erfcx_tail_q.cbegin(), x2);
    return 2.0 * pq.first / (x2 * pq.second);
  }
  const T u = 1.0 / x2;
  static constexpr std::array<double, 8> series_coefficients{
      1.0,     -1.5,       3.75,        -13.125,
      59.0625, -324.84375, 2111.484375, -15836.1328125};
  const T series = erfcx_paired_horner(series_coefficients, u, u * u);
  return -INV_SQRT_PI * u * series;
}

/**
 * Cody (1969) second-interval rational, valid for `0.46875 <= x <= 4`.
 *
 * Cody writes this interval as `erfc(x) = exp(-x*x) * R(x)`, so `erfcx`
 * is `R(x)` with the exponential cancelled analytically: no `exp` and no
 * `erfc` call.
 *
 * @param y argument, `0.46875 <= y <= 4`
 * @return scaled complementary error function
 */
inline double erfcx_cody_middle(double y) {
  const double y2 = y * y;
  return erfcx_paired_horner(erfcx_middle_p, y, y2)
         / erfcx_paired_horner(erfcx_middle_q, y, y2);
}

/**
 * Degree-18 polynomial for `|x| < 0.46875`.
 *
 * A Chebyshev-economized expansion of `erfcx` about zero, so it needs no
 * library call and no branch on the sign of `x`. Its low-order coefficients
 * reproduce the Maclaurin series of `erfcx` exactly (`1`, `-2/sqrt(pi)`,
 * `1`, ...), which is a useful check that the fit is right.
 *
 * Evaluated by Estrin's method: four independent degree-3 blocks are
 * combined in `x^4` and `x^8`.
 *
 * @param x argument, `|x| < 0.46875`
 * @return scaled complementary error function
 */
inline double erfcx_small(double x) {
  static constexpr std::array<double, 19> c{1.0,
                                            -1.12837916709551256,
                                            1.0,
                                            -7.52252778063651983e-01,
                                            4.99999999999992839e-01,
                                            -3.00901111227312890e-01,
                                            1.66666666667239644e-01,
                                            -8.59717459974174147e-02,
                                            4.16666666458337179e-02,
                                            -1.91048337772546720e-02,
                                            8.33333374332981443e-03,
                                            -3.47359067853470795e-03,
                                            1.38888415444527033e-03,
                                            -5.34506929034156810e-04,
                                            1.98445679338826757e-04,
                                            -7.08163358203131886e-05,
                                            2.46655529768908249e-05,
                                            -9.35890030086883823e-06,
                                            3.05977060678449757e-06};
  const double x2 = x * x;
  const double x4 = x2 * x2;
  const double x8 = x4 * x4;
  std::array<double, 4> quad;
  for (int i = 0; i < 4; ++i) {
    const int j = 4 * i + 2;
    quad[i] = (c[j] + c[j + 1] * x) + (c[j + 2] + c[j + 3] * x) * x2;
  }
  const double rest
      = (quad[0] + quad[1] * x4) + (quad[2] + quad[3] * x4 + c[18] * x8) * x8;
  return (c[0] + c[1] * x) + x2 * rest;
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
 * library function: a direct approximation of `erfcx` needs no `erfc`.
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
 *   into roughly `x * x * eps`. Below `x = -6.1`, `erfcx(-x)` is under
 *   `eps/2` of `2*exp(x*x)`, so the subtraction is skipped.
 *
 * W. J. Cody, Math. Comp. 23(107):631-638 (1969).
 * https://doi.org/10.1090/S0025-5718-1969-0247736-4
 *
 * The crossover at 4 is the same one R's `pnorm` uses: it switches to the
 * Cody tail form at `|y| > M_SQRT_32`, and since the argument here is the
 * `erfc` argument `y / sqrt(2)`, that is `sqrt(32) / sqrt(2) = 4` exactly.
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
