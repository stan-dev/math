// Original code derived from TensorFlow Probability and is distributed here
// under the Apache License 2.0 (licenses/tensorflow-probability-license.txt)
//    Copyright 2022 The TensorFlow Probability Authors.
//    tensorflow_probability/python/math/special.py, betainc partials
// Secondary code copyright by its author and is distributed here
// under the BSD-3 license (LICENSE.md)

#ifndef STAN_MATH_PRIM_FUN_GRAD_REG_INC_BETA_HPP
#define STAN_MATH_PRIM_FUN_GRAD_REG_INC_BETA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

namespace internal {

/**
 * Partial derivatives of the regularized incomplete beta function
 * I_x(a, b) with respect to a and b, by the power series
 *
 *   I_x(a, b) = x^a 2F1(a, 1 - b; a + 1; x) / (a B(a, b)),
 *
 * differentiated term by term. Used where at least one of these holds:
 *   x <  a / (a + b),  b x <= 1,        x <= 0.95
 *   x >= a / (a + b),  a (1 - x) <= 1,  x >= 0.05
 * For x >= a / (a + b) the symmetry relation I_x(a, b) = 1 - I_{1-x}(b, a)
 * is applied first. Stops when a term is below eps / a in magnitude, or
 * after 600 terms.
 *
 * @tparam T scalar type
 * @param a first shape
 * @param b second shape
 * @param x upper limit, in (0, 1)
 * @param digamma_a digamma(a)
 * @param digamma_b digamma(b)
 * @param digamma_ab digamma(a + b)
 * @param[out] g1 d/da I_x(a, b)
 * @param[out] g2 d/db I_x(a, b)
 */
template <typename T>
inline void inc_beta_partials_power_series(T a, T b, T x, T digamma_a,
                                           T digamma_b, T digamma_ab, T& g1,
                                           T& g2) {
  using std::exp;
  using std::log;
  const bool use_symmetry_relation
      = value_of_rec(x) >= value_of_rec(a) / value_of_rec(a + b);
  if (use_symmetry_relation) {
    const T a_orig = a;
    a = b;
    b = a_orig;
    x = 1.0 - x;
    const T digamma_a_orig = digamma_a;
    digamma_a = digamma_b;
    digamma_b = digamma_a_orig;
  }

  constexpr int max_iterations = 600;
  const double tolerance
      = std::numeric_limits<double>::epsilon() / value_of_rec(a);

  // series for 2F1(a, 1 - b; a + 1; x) / a and its a and b derivatives
  T product = 1.0;
  T series_sum = inv(a);
  T product_grad_b = 0.0;
  T da = -inv(square(a));
  T db = 0.0;
  for (int n = 1; n <= max_iterations; ++n) {
    const T x_div_n = x / n;
    const T factor = (n - b) * x_div_n;
    const T apn = a + n;
    const T new_product = product * factor;
    const T term = new_product / apn;
    product_grad_b = factor * product_grad_b - product * x_div_n;
    const T db_increment = product_grad_b / apn;
    da -= new_product / square(apn);
    db += db_increment;
    product = new_product;
    series_sum += term;
    // Differs from the TensorFlow Probability source: the stop test also
    // requires the d/db increment to be small. For a positive integer b the
    // value series is a polynomial whose term is exactly 0 at n = b, while
    // the d/db series continues.
    if (fabs(value_of_rec(term)) <= tolerance
        && fabs(value_of_rec(db_increment)) <= tolerance) {
      break;
    }
  }

  const T normalization = exp(a * log(x) - lbeta(a, b));
  const T grad_a
      = normalization * (da + series_sum * (digamma_ab - digamma_a + log(x)));
  const T grad_b = normalization * (db + series_sum * (digamma_ab - digamma_b));

  if (use_symmetry_relation) {
    g1 = -grad_b;
    g2 = -grad_a;
  } else {
    g1 = grad_a;
    g2 = grad_b;
  }
}

/**
 * One step of the modified Lentz method for the continued fraction of the
 * incomplete beta function, carrying the derivatives of the ratios with
 * respect to a and b. Thompson and Barnett (1986), appendix.
 *
 * @tparam T scalar type
 * @param d partial numerator d_n
 * @param dd_a d/da of d_n
 * @param dd_b d/db of d_n
 * @param[in,out] C ratio of successive numerators A_n / A_{n-1}
 * @param[in,out] D ratio of successive denominators B_{n-1} / B_n
 * @param[in,out] h convergent A_n / B_n
 * @param[in,out] dC_a d/da of C
 * @param[in,out] dC_b d/db of C
 * @param[in,out] dD_a d/da of D
 * @param[in,out] dD_b d/db of D
 * @param[in,out] dh_a d/da of h
 * @param[in,out] dh_b d/db of h
 * @return the multiplicative update of h at this step
 */
template <typename T>
inline T inc_beta_lentz_step(const T& d, const T& dd_a, const T& dd_b, T& C,
                             T& D, T& h, T& dC_a, T& dC_b, T& dD_a, T& dD_b,
                             T& dh_a, T& dh_b) {
  const double small = std::sqrt(std::numeric_limits<double>::min());

  T C_new = 1.0 + d / C;
  if (fabs(value_of_rec(C_new)) < small) {
    C_new = small;
  }
  T D_new = 1.0 + d * D;
  if (fabs(value_of_rec(D_new)) < small) {
    D_new = small;
  }
  D_new = inv(D_new);
  const T delta = C_new * D_new;
  const T h_new = h * delta;

  const T dC_a_new = (dd_a * C - d * dC_a) / square(C);
  const T dC_b_new = (dd_b * C - d * dC_b) / square(C);
  const T dD_a_new = -(dd_a * D + d * dD_a) * square(D_new);
  const T dD_b_new = -(dd_b * D + d * dD_b) * square(D_new);
  const T dh_a_new = dh_a * delta + h * dC_a_new * D_new + h * dD_a_new * C_new;
  const T dh_b_new = dh_b * delta + h * dC_b_new * D_new + h * dD_b_new * C_new;

  C = C_new;
  D = D_new;
  h = h_new;
  dC_a = dC_a_new;
  dC_b = dC_b_new;
  dD_a = dD_a_new;
  dD_b = dD_b_new;
  dh_a = dh_a_new;
  dh_b = dh_b_new;
  return delta;
}

/**
 * Partial derivatives of the regularized incomplete beta function
 * I_x(a, b) with respect to a and b, by the continued fraction
 * https://dlmf.nist.gov/8.17.E22 evaluated with the modified Lentz method,
 * with the partial numerators https://dlmf.nist.gov/8.17.E23 differentiated
 * in a and b and the derivatives carried through the recurrence. Used
 * where the power series region does not apply. The continued fraction
 * converges rapidly for x < (a - 1) / (a + b - 2); otherwise the symmetry
 * relation I_x(a, b) = 1 - I_{1-x}(b, a) is applied first. Stops when the
 * multiplicative update is within 3 eps of 1, or after 300 double steps.
 *
 * @tparam T scalar type
 * @param a first shape
 * @param b second shape
 * @param x upper limit, in (0, 1)
 * @param digamma_a digamma(a)
 * @param digamma_b digamma(b)
 * @param digamma_ab digamma(a + b)
 * @param[out] g1 d/da I_x(a, b)
 * @param[out] g2 d/db I_x(a, b)
 */
template <typename T>
inline void inc_beta_partials_continued_fraction(T a, T b, T x, T digamma_a,
                                                 T digamma_b, T digamma_ab,
                                                 T& g1, T& g2) {
  using std::exp;
  using std::log;
  const double a_d = value_of_rec(a);
  const double b_d = value_of_rec(b);
  // IEEE division: (a - 1) / 0 is +-inf and the comparison is well defined
  const bool use_symmetry_relation
      = value_of_rec(x) >= (a_d - 1.0) / (a_d + b_d - 2.0);
  if (use_symmetry_relation) {
    const T a_orig = a;
    a = b;
    b = a_orig;
    x = 1.0 - x;
    const T digamma_a_orig = digamma_a;
    digamma_a = digamma_b;
    digamma_b = digamma_a_orig;
  }

  constexpr int max_iterations = 300;
  const double tolerance = 3.0 * std::numeric_limits<double>::epsilon();
  const double small = std::sqrt(std::numeric_limits<double>::min());

  const T apb = a + b;
  const T ap1 = a + 1.0;

  // initialization and first step of the modified Lentz method
  T C = 1.0;
  T D = 1.0 - apb * x / ap1;
  if (fabs(value_of_rec(D)) < small) {
    D = small;
  }
  D = inv(D);
  T h = D;
  const T dD_denom = square(x * apb - ap1);
  T dD_a = (1.0 - b) * x / dD_denom;
  T dD_b = ap1 * x / dD_denom;
  T dC_a = 0.0;
  T dC_b = 0.0;
  T dh_a = dD_a;
  T dh_b = dD_b;

  for (int m = 1; m <= max_iterations; ++m) {
    const double dh_a_prev = value_of_rec(dh_a);
    const double dh_b_prev = value_of_rec(dh_b);
    // even partial numerator d_{2m}, https://dlmf.nist.gov/8.17.E23
    {
      const T a_plus_2m = a + 2.0 * m;
      const T a_plus_2m_minus_one = a_plus_2m - 1.0;
      const T denominator = a_plus_2m * a_plus_2m_minus_one;
      const T dd_b = m * x / denominator;
      const T d = dd_b * (b - m);
      const T dd_a = -d * (a_plus_2m + a_plus_2m_minus_one) / denominator;
      inc_beta_lentz_step(d, dd_a, dd_b, C, D, h, dC_a, dC_b, dD_a, dD_b, dh_a,
                          dh_b);
    }
    // odd partial numerator d_{2m+1}
    T delta;
    {
      const T a_plus_m = a + m;
      const T a_plus_2m = a_plus_m + m;
      const T a_plus_2m_plus_one = a_plus_2m + 1.0;
      const T a_plus_b_plus_m = a_plus_m + b;
      const T denominator = a_plus_2m * a_plus_2m_plus_one;
      const T dd_b = -a_plus_m * x / denominator;
      const T d = dd_b * a_plus_b_plus_m;
      const T dd_a = -d * ((a_plus_2m + a_plus_2m_plus_one) / denominator)
                     - x * (2.0 * a_plus_m + b) / denominator;
      delta = inc_beta_lentz_step(d, dd_a, dd_b, C, D, h, dC_a, dC_b, dD_a,
                                  dD_b, dh_a, dh_b);
    }
    // Differs from the TensorFlow Probability source: the stop test also
    // requires the two derivatives to have converged. For a positive
    // integer b the partial numerator d_{2m} is exactly 0 at m = b; the
    // convergent is then exact and delta is 1, but the derivatives still
    // need the tail of the fraction.
    const double dh_a_now = value_of_rec(dh_a);
    const double dh_b_now = value_of_rec(dh_b);
    if (fabs(value_of_rec(delta) - 1.0) < tolerance
        && fabs(dh_a_now - dh_a_prev) <= tolerance * fabs(dh_a_now)
        && fabs(dh_b_now - dh_b_prev) <= tolerance * fabs(dh_b_now)) {
      break;
    }
  }

  const T normalization
      = exp(a * log(x) + b * log1p(-x) - log(a) - lbeta(a, b));
  const T grad_a
      = normalization * (dh_a + h * (log(x) - inv(a) + digamma_ab - digamma_a));
  const T grad_b
      = normalization * (dh_b + h * (log1p(-x) + digamma_ab - digamma_b));

  if (use_symmetry_relation) {
    g1 = -grad_b;
    g2 = -grad_a;
  } else {
    g1 = grad_a;
    g2 = grad_b;
  }
}

}  // namespace internal

/**
 * Computes the gradients of the regularized incomplete beta
 * function.  Specifically, this function computes gradients of
 * <code>inc_beta(a, b, z)</code>, with respect to the arguments
 * <code>a</code> and <code>b</code>.
 *
 * The algorithm is the one of Boik and Robinson-Cox (1998), "Derivatives
 * of the Incomplete Beta Function", Journal of Statistical Software 3(1),
 * in the form implemented by TensorFlow Probability: a power series near
 * the endpoints and a continued fraction elsewhere, both differentiated
 * term by term, with the prefactor evaluated in log space and the
 * symmetry relation I_z(a, b) = 1 - I_{1-z}(b, a) applied where the
 * expansion converges faster on the other side.
 *
 * The gradients are 0 at z = 0 and z = 1.
 *
 * @tparam T type of arguments
 * @param[out] g1 partial derivative of <code>inc_beta(a, b, z)</code>
 * with respect to <code>a</code>
 * @param[out] g2 partial derivative of <code>inc_beta(a, b,
 * z)</code> with respect to <code>b</code>
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @param[in] digammaA the value of <code>digamma(a)</code>
 * @param[in] digammaB the value of <code>digamma(b)</code>
 * @param[in] digammaSum the value of <code>digamma(a + b)</code>
 * @param[in] betaAB the value of <code>beta(a, b)</code>; accepted for
 * interface compatibility and not used, the normalization is evaluated in
 * log space
 */
template <typename T>
inline void grad_reg_inc_beta(T& g1, T& g2, const T& a, const T& b, const T& z,
                              const T& digammaA, const T& digammaB,
                              const T& digammaSum, const T& betaAB) {
  const double a_d = value_of_rec(a);
  const double b_d = value_of_rec(b);
  const double z_d = value_of_rec(z);

  if (z_d == 0.0 || z_d == 1.0) {
    g1 = 0.0;
    g2 = 0.0;
    return;
  }

  const double mean = a_d / (a_d + b_d);
  const bool use_power_series
      = (z_d < mean && b_d * z_d <= 1.0 && z_d <= 0.95)
        || (z_d >= mean && a_d * (1.0 - z_d) <= 1.0 && z_d >= 0.05);
  if (use_power_series) {
    internal::inc_beta_partials_power_series(a, b, z, digammaA, digammaB,
                                             digammaSum, g1, g2);
  } else {
    internal::inc_beta_partials_continued_fraction(a, b, z, digammaA, digammaB,
                                                   digammaSum, g1, g2);
  }
}

}  // namespace math
}  // namespace stan
#endif
