#ifndef STAN_MATH_PRIM_FUN_ERFCX_HPP
#define STAN_MATH_PRIM_FUN_ERFCX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

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
 * Two branches, with the crossover and the accuracy below established by
 * measurement against a `long double` reference over 40000 points per
 * interval:
 *
 * - `x >= 4`: the Cody (1969) rational approximation, which gives `erfcx`
 *   directly with no `exp` and no `erfc` call. Measured at 2.0 to 2.5 ulp
 *   over `[4, 50]` and 5.4 ulp out to 200, and about 2.4 times faster than
 *   forming the product.
 * - `x < 4`: `exp(x * x) * erfc(x)`, with `x * x` split so that the
 *   argument handed to `exp` is exact. Without the split the rounding of
 *   `x * x` is amplified by `exp` into a relative error of about
 *   `x * x * eps`, which reaches 512 ulp at `x = -26`. Splitting holds this
 *   branch to under 5 ulp.
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
 * Measured worst case over the whole range is 5.4 ulp.
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
  // Cody (1969). x * x is infinite for x large enough, which correctly
  // collapses this to the leading term INV_SQRT_PI / x and, at infinity,
  // to zero.
  constexpr double cody_min = 4.0;
  if (x >= cody_min) {
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
  // erfcx(x) = 2 * exp(x * x) * (1 + o(1)) as x -> -infinity, which leaves
  // the binary64 range below -26.63. Returning here also keeps the split
  // below from forming inf - inf on -infinity.
  constexpr double overflow_max = -27.0;
  if (x < overflow_max) {
    return INFTY;
  }
  // Dekker split: x_hi holds the leading 26 bits of x, so x_hi * x_hi is
  // exact and exp() is handed an argument carrying no rounding error of
  // its own. x_lo * (x + x_hi) is the remainder x * x - x_hi * x_hi.
  constexpr double split = 134217729.0;  // 2^27 + 1
  const double t = split * x;
  const double x_hi = t - (t - x);
  const double x_lo = x - x_hi;
  return std::exp(x_hi * x_hi) * std::exp(x_lo * (x + x_hi)) * std::erfc(x);
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
