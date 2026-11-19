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

/**
 * Above this the tail is the leading term 1 / (x sqrt(pi)) alone. The
 * correction is 1 / (2 x^2) of that term, which at 2^27 = 134217728 is
 * exactly eps / 8, so from here on the single division is the more accurate
 * result. The rational's x^13 denominator would overflow above about 2^79,
 * so some such cut is also what keeps erfcx finite out to +infinity.
 */
inline constexpr double erfcx_tail_leading_only = 134217728.0;

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
 * Horner evaluation of the two Cody tail polynomials `erfcx_tail_p` and
 * `erfcx_tail_q` at `t`.
 *
 * The arrays are stored in ascending powers of `u = 1 / x^2`. Read front to
 * back that is descending powers of `s = x^2`, so the default walk is Horner
 * in `s`. With `Reverse` the arrays are read back to front, which is Horner
 * in `u`.
 *
 * @tparam Reverse read the coefficients back to front
 * @tparam T scalar type
 * @param t argument, `x^2` or with `Reverse` `1 / x^2`
 * @return numerator and denominator
 */
template <bool Reverse = false, typename T>
inline constexpr std::pair<T, T> erfcx_tail_polynomials(const T& t) {
  constexpr std::size_t N = erfcx_tail_p.size();
  constexpr auto at = [](std::size_t i) { return Reverse ? N - 1 - i : i; };
  T numerator = erfcx_tail_p[at(0)];
  T denominator = erfcx_tail_q[at(0)];
  for (std::size_t i = 1; i < N; ++i) {
    numerator = numerator * t + erfcx_tail_p[at(i)];
    denominator = denominator * t + erfcx_tail_q[at(i)];
  }
  return {numerator, denominator};
}

/**
 * Numerator and denominator of the Cody (1969) third-interval rational,
 * valid for `x >= 4`.
 *
 * `erfcx(x) = (INV_SQRT_PI + u * C(u)) / x` with `u = 1 / x^2` and `C` the
 * Cody correction. `u * C(u)` is formed from the coefficients in `x^2`, so
 * `1 / x^2` is never formed, and the sum is put over the common denominator
 * `x^3 Q(x^2)` so the quotient takes one division. Below
 * `erfcx_tail_leading_only` that denominator is at most about `2^351`.
 *
 * @tparam T scalar type
 * @param x argument, `x >= 4`
 * @return numerator and denominator of `erfcx(x)`
 */
template <typename T>
inline constexpr std::pair<T, T> erfcx_cody_tail_fraction(const T& x) {
  if (x > erfcx_tail_leading_only) {
    return {T(INV_SQRT_PI), x};
  }
  const T s = x * x;
  const auto pq = erfcx_tail_polynomials(s);
  const T s_q = s * pq.second;
  return {INV_SQRT_PI * s_q + pq.first, x * s_q};
}

/**
 * Cody (1969) third-interval rational, valid for `x >= 4`.
 *
 * @tparam T scalar type
 * @param x argument, `x >= 4`
 * @return scaled complementary error function
 */
template <typename T>
inline constexpr T erfcx_cody_tail(const T& x) {
  const auto pq = erfcx_cody_tail_fraction(x);
  return pq.first / pq.second;
}

static constexpr std::array<double, 8> erfcx_deriv_c{
    1.0, -1.5, 3.75, -13.125, 59.0625, -324.84375, 2111.484375, -15836.1328125};
// Horner in u^2 over the pairs (c[6], c[7]), (c[4], c[5]), ..., (c[0], c[1])

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
inline constexpr T erfcx_derivative(const T& x, const T& value) {
  if (x < 4.0) {
    return 2.0 * x * value - TWO_OVER_SQRT_PI;
  }
  const T x2 = x * x;
  if (x < 30.0) {
    const auto pq = erfcx_tail_polynomials(x2);
    return 2.0 * pq.first / (x2 * pq.second);
  }
  const T u = 1.0 / x2;
  const T u2 = u * u;
  constexpr auto&& c = erfcx_deriv_c;
  T series = c[6] + c[7] * u;
  for (int j = 4; j >= 0; j -= 2) {
    series = series * u2 + (c[j] + c[j + 1] * u);
  }
  return -INV_SQRT_PI * u * series;
}

/**
 * Numerator and denominator of the Cody (1969) second-interval rational,
 * valid for `0.46875 <= x <= 4`.
 *
 * Cody writes this interval as `erfc(x) = exp(-x*x) * R(x)`, so `erfcx`
 * is `R(x)` with the exponential cancelled analytically: no `exp` and no
 * `erfc` call.
 *
 * Both degree-8 polynomials are Horner in `y^2` over adjacent coefficient
 * pairs, seeded by the lone top coefficient, so the `c[j + 1] * y` products
 * sit off the dependency chain.
 *
 * @tparam T scalar type
 * @param y argument, `0.46875 <= y <= 4`
 * @return numerator and denominator of `erfcx(y)`
 */
template <typename T>
inline constexpr std::pair<T, T> erfcx_cody_middle_fraction(const T& y) {
  const T y2 = y * y;
  T p = erfcx_middle_p[8];
  T q = erfcx_middle_q[8];
  for (int j = 6; j >= 0; j -= 2) {
    p = p * y2 + (erfcx_middle_p[j] + erfcx_middle_p[j + 1] * y);
    q = q * y2 + (erfcx_middle_q[j] + erfcx_middle_q[j + 1] * y);
  }
  return {p, q};
}

/**
 * Cody (1969) second-interval rational, valid for `0.46875 <= x <= 4`.
 *
 * @tparam T scalar type
 * @param y argument, `0.46875 <= y <= 4`
 * @return scaled complementary error function
 */
template <typename T>
inline constexpr T erfcx_cody_middle(const T& y) {
  const auto pq = erfcx_cody_middle_fraction(y);
  return pq.first / pq.second;
}

// R's coefficients by the power of x they multiply
static constexpr std::array<double, 4> small_r0_3{1.0, -7.52252778063651983e-01,
                                                  4.99999999999992839e-01,
                                                  -3.00901111227312890e-01};
static constexpr std::array<double, 4> small_r4_7{
    1.66666666667239644e-01, -8.59717459974174147e-02, 4.16666666458337179e-02,
    -1.91048337772546720e-02};
static constexpr std::array<double, 4> small_r8_11{
    8.33333374332981443e-03, -3.47359067853470795e-03, 1.38888415444527033e-03,
    -5.34506929034156810e-04};
static constexpr std::array<double, 4> small_r12_15{
    1.98445679338826757e-04, -7.08163358203131886e-05, 2.46655529768908249e-05,
    -9.35890030086883823e-06};
static constexpr double small_r16 = 3.05977060678449757e-06;

/**
 * Degree-18 polynomial for `|x| < 0.46875`.
 *
 * A Chebyshev-economized expansion of `erfcx` about zero, so it needs no
 * library call and no branch on the sign of `x`. Its low-order coefficients
 * reproduce the Maclaurin series of `erfcx` exactly (`1`, `-2/sqrt(pi)`,
 * `1`, ...), which is a useful check that the fit is right.
 *
 * Written as `erfcx(x) = 1 - 2 x / sqrt(pi) + x^2 R(x)`: the first two terms
 * are the exact Maclaurin coefficients and `R` is the degree-16 remainder of
 * the fit. Adding `x^2 R` to the exact pair last is the accurate order for
 * small `x`. `R` is evaluated by Estrin's method: four cubic blocks in `x`,
 * combined in `x^4` and `x^8`, with the degree-16 coefficient on top.
 *
 * @tparam T scalar type
 * @param x argument, `|x| < 0.46875`
 * @return scaled complementary error function
 */
template <typename T>
inline constexpr T erfcx_small(const T& x) {
  constexpr double c0 = 1.0;
  constexpr double c1 = -1.12837916709551256;
  constexpr auto&& r0_3 = small_r0_3;
  constexpr auto&& r4_7 = small_r4_7;
  constexpr auto&& r8_11 = small_r8_11;
  constexpr auto&& r12_15 = small_r12_15;
  constexpr auto&& r16 = small_r16;
  const T x2 = x * x;
  const T x4 = x2 * x2;
  const T x8 = x4 * x4;
  const T q0 = (r0_3[0] + r0_3[1] * x) + (r0_3[2] + r0_3[3] * x) * x2;
  const T q1 = (r4_7[0] + r4_7[1] * x) + (r4_7[2] + r4_7[3] * x) * x2;
  const T q2 = (r8_11[0] + r8_11[1] * x) + (r8_11[2] + r8_11[3] * x) * x2;
  const T q3 = (r12_15[0] + r12_15[1] * x) + (r12_15[2] + r12_15[3] * x) * x2;
  const T rest = (q0 + q1 * x4) + (q2 + q3 * x4 + r16 * x8) * x8;
  return (c0 + c1 * x) + x2 * rest;
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
  } else if (-x >= cody_min) {
    return two_exp_x2 - internal::erfcx_cody_tail(-x);
  } else {
    return two_exp_x2 - internal::erfcx_cody_middle(-x);
  }
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
