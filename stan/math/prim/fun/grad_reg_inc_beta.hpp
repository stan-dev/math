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
#include <cstdint>
#include <limits>
#include <utility>

namespace stan {
namespace math {

/**
 * Selects which partial derivatives of the regularized incomplete beta
 * function <code>grad_reg_inc_beta</code> computes: <code>G1</code> for the
 * derivative with respect to a, <code>G2</code> for the derivative with
 * respect to b, <code>BOTH</code> for both.
 */
enum class GRAD_RIB : uint8_t { G1 = 0, G2 = 1, BOTH = 2 };

namespace internal {

/**
 * Partial derivatives of the regularized incomplete beta function I_x(a, b)
 * with respect to a and b, computed from its power series.
 *
 * The series is
 *
 *   I_x(a, b) = N S,  N = x^a / B(a, b),  S = sum_{n >= 0} t_n,
 *   t_n = P_n / (a + n),  P_0 = 1,  P_n = P_{n-1} (n - b) x / n,
 *
 * where S is 2F1(a, 1 - b; a + 1; x) / a. The derivatives of each term are
 *
 *   dt_n/da = -t_n / (a + n),  dt_n/db = (dP_n/db) / (a + n),
 *
 * and dP_n/db follows from the recurrence for P_n by the product rule. With
 * log N = a log x - log B(a, b) and psi the digamma function,
 *
 *   dI/da = N (dS/da + S (log x - psi(a) + psi(a + b))),
 *   dI/db = N (dS/db + S (psi(a + b) - psi(b))).
 *
 * Only the sums needed for the requested partials are accumulated.
 *
 * The loop stops once a term is below eps / a in magnitude, or after 600
 * terms. For a positive integer b, P_n is exactly 0 for n >= b but dP_n/db is
 * not, so the sum for dS/db keeps changing after S has converged. The stop
 * test therefore also requires dt_n/db to be small whenever that sum
 * matters. That is the case when d/db is requested, and also when b is an
 * autodiff type, because autodiff then differentiates the d/da result in b
 * and needs the same terms.
 *
 * @tparam grads_requested which partials to compute
 * @tparam flip_sign true when the caller evaluates its I_z(a, b) through the
 *   symmetry relation I_z(a, b) = 1 - I_{1-z}(b, a) and so passes
 *   (b, a, 1 - z). The results are then negated, and for <code>BOTH</code>
 *   they are returned in the caller's order of a and b.
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of x
 * @tparam T4 type of digamma_a
 * @tparam T5 type of digamma_b
 * @tparam T6 type of digamma_ab
 * @param a first shape
 * @param b second shape
 * @param x upper limit, in (0, 1)
 * @param digamma_a digamma(a)
 * @param digamma_b digamma(b)
 * @param digamma_ab digamma(a + b)
 * @return d/da I_x(a, b) for <code>G1</code>, d/db for <code>G2</code>, and
 *   the pair (d/da, d/db) for <code>BOTH</code>. With <code>flip_sign</code>
 *   these are negated and the pair is returned as (-d/db, -d/da), which are
 *   the caller's (d/da, d/db).
 */
template <GRAD_RIB grads_requested, bool flip_sign, typename T1, typename T2,
          typename T3, typename T4, typename T5, typename T6>
inline auto inc_beta_power_series_kernel(const T1& a, const T2& b, const T3& x,
                                         const T4& digamma_a,
                                         const T5& digamma_b,
                                         const T6& digamma_ab) {
  using ret_t = return_type_t<T1, T2, T3, T4, T5, T6>;
  constexpr int max_iterations = 600;
  const double tolerance
      = std::numeric_limits<double>::epsilon() / value_of_rec(a);

  // series for 2F1(a, 1 - b; a + 1; x) / a and its a and b derivatives
  ret_t product = 1.0;
  ret_t series_sum = inv(a);
  ret_t da = -series_sum / a;
  ret_t db = 0.0;
  ret_t product_grad_b = 0.0;
  for (int n = 1; n <= max_iterations; ++n) {
    const auto x_div_n = x / n;
    const auto factor = (n - b) * x_div_n;
    const auto apn = a + n;
    const auto new_product = product * factor;
    const auto term = new_product / apn;
    bool db_converged = true;
    if constexpr (grads_requested != GRAD_RIB::G1 || is_autodiff_v<T2>) {
      product_grad_b = factor * product_grad_b - product * x_div_n;
      const auto db_increment = product_grad_b / apn;
      if constexpr (grads_requested != GRAD_RIB::G1) {
        db += db_increment;
      }
      db_converged = fabs(value_of_rec(db_increment)) <= tolerance;
    }
    if constexpr (grads_requested != GRAD_RIB::G2) {
      da -= term / apn;
    }
    product = new_product;
    series_sum += term;
    if (fabs(value_of_rec(term)) <= tolerance && db_converged) {
      break;
    }
  }

  // the symmetry relation negates the partials
  constexpr double sign = flip_sign ? -1.0 : 1.0;
  const auto log_x = log(x);
  const auto normalize = sign * exp(a * log_x - lbeta(a, b));
  if constexpr (grads_requested == GRAD_RIB::G1) {
    return normalize * (da + series_sum * (digamma_ab - digamma_a + log_x));
  } else if constexpr (grads_requested == GRAD_RIB::G2) {
    return normalize * (db + series_sum * (digamma_ab - digamma_b));
  } else {
    const auto grad_a
        = normalize * (da + series_sum * (digamma_ab - digamma_a + log_x));
    const auto grad_b
        = normalize * (db + series_sum * (digamma_ab - digamma_b));
    if constexpr (flip_sign) {
      // the caller swapped a and b, so its (d/da, d/db) are (grad_b, grad_a)
      return std::pair{grad_b, grad_a};
    } else {
      return std::pair{grad_a, grad_b};
    }
  }
}

/**
 * Partial derivatives of the regularized incomplete beta function I_x(a, b)
 * with respect to a and b, computed from its continued fraction.
 *
 * The continued fraction https://dlmf.nist.gov/8.17.E22 is
 *
 *   I_x(a, b) = Q h,  Q = x^a (1 - x)^b / (a B(a, b)),
 *   h = 1 / (1 + d_1 / (1 + d_2 / (1 + ...))),
 *
 * with the partial numerators https://dlmf.nist.gov/8.17.E23
 *
 *   d_1 = -(a + b) x / (a + 1),
 *   d_{2m} = m (b - m) x / ((a + 2m - 1)(a + 2m)),
 *   d_{2m+1} = -(a + m)(a + b + m) x / ((a + 2m)(a + 2m + 1)).
 *
 * Each d_n is num / den with den = (a + n - 1)(a + n), so its a derivative is
 * (num_a - d_n den_a) / den, with den_a = 2 (a + n) - 1. The even numerators
 * do not depend on a, and the odd ones have num_a = -(2 (a + m) + b) x. Every
 * numerator is linear in b, so the b derivative is num_b / den.
 *
 * h is evaluated with the modified Lentz method. With h_n = A_n / B_n the
 * n-th convergent, it keeps the ratios C = A_n / A_{n-1} and
 * D = B_{n-1} / B_n, and step n updates them as
 *
 *   C <- 1 + d_n / C,  D <- 1 / (1 + d_n D),  h <- h C D.
 *
 * The first step holds C at 1, so D = 1 / (1 + d_1) and h = D. When
 * 1 + d_n / C or 1 + d_n D is smaller in magnitude than sqrt(DBL_MIN), it is
 * replaced by that value, so that nothing divides by 0.
 *
 * The derivatives are carried through the same steps, following Thompson and
 * Barnett (1986), appendix. They are kept as logarithmic derivatives
 * c = (dC/da) / C, e = (dD/da) / D, and r = (dh/da) / h, and likewise for b.
 * Differentiating the updates gives
 *
 *   c <- (dd_n/da - d_n c) / (C C_new),
 *   e <- -(dd_n/da + d_n e) D D_new,
 *   r <- r + c + e,
 *
 * where C, D, c, and e on the right are the values before the step and
 * C_new, D_new are the values after it. At the end, with psi the digamma
 * function,
 *
 *   dI/da = I (r_a + log x - 1 / a + psi(a + b) - psi(a)),
 *   dI/db = I (r_b + log(1 - x) + psi(a + b) - psi(b)).
 *
 * Only the requested derivatives are carried.
 *
 * The test for stopping runs after each odd n. The loop stops once C D is
 * within 3 eps of 1 and each carried derivative dh = h r has changed by at
 * most a relative 3 eps since the previous test, or after n = 601. For
 * integer b the partial numerator d_{2b} is 0, so from that step on h and its
 * a derivative no longer change. The b derivative of d_{2b} is b x / den,
 * which is not 0, so dh/db is still changing. The b derivative is therefore
 * tested whenever it matters, which is when d/db is requested and also when b
 * is an autodiff type, as in the power series.
 *
 * @tparam grads_requested which partials to compute
 * @tparam flip_sign true when the caller evaluates its I_z(a, b) through the
 *   symmetry relation I_z(a, b) = 1 - I_{1-z}(b, a) and so passes
 *   (b, a, 1 - z). The results are then negated, and for <code>BOTH</code>
 *   they are returned in the caller's order of a and b.
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of x
 * @tparam T4 type of digamma_a
 * @tparam T5 type of digamma_b
 * @tparam T6 type of digamma_ab
 * @param a first shape
 * @param b second shape
 * @param x upper limit, in (0, 1)
 * @param digamma_a digamma(a)
 * @param digamma_b digamma(b)
 * @param digamma_ab digamma(a + b)
 * @return d/da I_x(a, b) for <code>G1</code>, d/db for <code>G2</code>, and
 *   the pair (d/da, d/db) for <code>BOTH</code>. With <code>flip_sign</code>
 *   these are negated and the pair is returned as (-d/db, -d/da), which are
 *   the caller's (d/da, d/db).
 */
template <GRAD_RIB grads_requested, bool flip_sign, typename T1, typename T2,
          typename T3, typename T4, typename T5, typename T6>
inline auto inc_beta_continued_fraction_kernel(const T1& a, const T2& b,
                                               const T3& x, const T4& digamma_a,
                                               const T5& digamma_b,
                                               const T6& digamma_ab) {
  using ret_t = return_type_t<T1, T2, T3, T4, T5, T6>;
  constexpr int max_iterations = 300;
  const double tolerance = 3.0 * std::numeric_limits<double>::epsilon();
  const double small = std::sqrt(std::numeric_limits<double>::min());

  // first step, with C held at 1
  const auto ap1 = a + 1.0;
  const auto d_1 = -(a + b) * x / ap1;
  ret_t C = 1.0;
  ret_t D = 1.0 + d_1;
  if (fabs(value_of_rec(D)) < small) {
    D = small;
  }
  D = inv(D);
  ret_t h = D;
  ret_t c_a = 0.0;
  ret_t e_a = 0.0;
  ret_t r_a = 0.0;
  ret_t c_b = 0.0;
  ret_t e_b = 0.0;
  ret_t r_b = 0.0;
  if constexpr (grads_requested != GRAD_RIB::G2) {
    // d d_1 / da = -(x + d_1) / (a + 1)
    e_a = (x + d_1) / ap1 * D;
    r_a = e_a;
  }
  if constexpr (grads_requested != GRAD_RIB::G1 || is_autodiff_v<T2>) {
    // d d_1 / db = -x / (a + 1)
    e_b = x / ap1 * D;
    r_b = e_b;
  }

  // the stop test compares dh = h r after each odd step
  double dh_a_prev = value_of_rec(h) * value_of_rec(r_a);
  double dh_b_prev = value_of_rec(h) * value_of_rec(r_b);
  for (int n = 2; n <= 2 * max_iterations + 1; ++n) {
    // the partial numerator d_n and its b derivative dd_b
    const int m = n / 2;
    const bool odd = n % 2 == 1;
    const auto a_plus_n = a + n;
    const auto denominator = a_plus_n * (a_plus_n - 1.0);
    ret_t dd_b;
    ret_t d;
    if (odd) {
      dd_b = -(a + m) * x / denominator;
      d = dd_b * (a + b + m);
    } else {
      dd_b = m * x / denominator;
      d = dd_b * (b - m);
    }
    // one step of the modified Lentz method
    auto C_new = 1.0 + d / C;
    // C C_new = C + d, which does not wait on the division for C_new
    ret_t C_C_new = C + d;
    if (fabs(value_of_rec(C_new)) < small) {
      C_new = small;
      C_C_new = C * small;
    }
    auto D_new = 1.0 + d * D;
    if (fabs(value_of_rec(D_new)) < small) {
      D_new = small;
    }
    D_new = inv(D_new);
    // c_new = (dd - d c) / (C C_new) and e_new = -(dd + d e) D D_new
    const auto D_D_new = D * D_new;
    if constexpr (grads_requested != GRAD_RIB::G2) {
      auto dd_a = -d * (2.0 * a_plus_n - 1.0) / denominator;
      if (odd) {
        dd_a -= x * (2.0 * (a + m) + b) / denominator;
      }
      c_a = (dd_a - d * c_a) / C_C_new;
      e_a = -(dd_a + d * e_a) * D_D_new;
      r_a += c_a + e_a;
    }
    if constexpr (grads_requested != GRAD_RIB::G1 || is_autodiff_v<T2>) {
      c_b = (dd_b - d * c_b) / C_C_new;
      e_b = -(dd_b + d * e_b) * D_D_new;
      r_b += c_b + e_b;
    }
    C = C_new;
    D = D_new;
    const auto delta = C_new * D_new;
    h *= delta;
    if (odd) {
      const double h_d = value_of_rec(h);
      const double dh_a_now = h_d * value_of_rec(r_a);
      const double dh_b_now = h_d * value_of_rec(r_b);
      if (fabs(value_of_rec(delta) - 1.0) < tolerance
          && (grads_requested == GRAD_RIB::G2
              || fabs(dh_a_now - dh_a_prev) <= tolerance * fabs(dh_a_now))
          && ((grads_requested == GRAD_RIB::G1 && !is_autodiff_v<T2>)
              || fabs(dh_b_now - dh_b_prev) <= tolerance * fabs(dh_b_now))) {
        break;
      }
      dh_a_prev = dh_a_now;
      dh_b_prev = dh_b_now;
    }
  }

  // the symmetry relation negates the partials
  constexpr double sign = flip_sign ? -1.0 : 1.0;
  const auto log_x = log(x);
  const auto log1m_x = log1p(-x);
  // h times the normalization is I_x(a, b), shared by both partials
  const auto prefactor
      = sign * h * exp(a * log_x + b * log1m_x - log(a) - lbeta(a, b));
  if constexpr (grads_requested == GRAD_RIB::G1) {
    return ret_t(prefactor * (r_a + log_x - inv(a) + digamma_ab - digamma_a));
  } else if constexpr (grads_requested == GRAD_RIB::G2) {
    return ret_t(prefactor * (r_b + log1m_x + digamma_ab - digamma_b));
  } else {
    const ret_t grad_a
        = prefactor * (r_a + log_x - inv(a) + digamma_ab - digamma_a);
    const ret_t grad_b = prefactor * (r_b + log1m_x + digamma_ab - digamma_b);
    if constexpr (flip_sign) {
      // the caller swapped a and b, so its (d/da, d/db) are (grad_b, grad_a)
      return std::pair<ret_t, ret_t>{grad_b, grad_a};
    } else {
      return std::pair<ret_t, ret_t>{grad_a, grad_b};
    }
  }
}

}  // namespace internal

/**
 * Computes the requested partial derivatives of the regularized incomplete
 * beta function I_z(a, b) = <code>inc_beta(a, b, z)</code> with respect to
 * <code>a</code> and <code>b</code>.
 *
 * The method is that of Boik and Robinson-Cox (1998), "Derivatives of the
 * Incomplete Beta Function", Journal of Statistical Software 3(1), as
 * implemented in TensorFlow Probability. I_z(a, b) is expanded either as a
 * power series or as a continued fraction, and the expansion is
 * differentiated term by term. Only the requested derivatives are carried
 * through it. The normalizing factor is evaluated in log space, so it stays
 * finite where B(a, b) underflows.
 *
 * Either expansion can instead be applied to the right side of the symmetry
 * relation I_z(a, b) = 1 - I_{1-z}(b, a), which converges faster for some
 * arguments. The expansion is then evaluated for I_{1-z}(b, a), and the
 * partials of I_z(a, b) are its partials negated, with the roles of a and b
 * swapped.
 *
 * With mean = a / (a + b), the branch is chosen as follows.
 *
 * - The power series for I_z(a, b) if z < mean, b z <= 1, and z <= 0.5.
 * - The power series for I_{1-z}(b, a) if z >= mean, a (1 - z) <= 1, and
 *   z >= 0.5.
 * - Otherwise, the continued fraction for I_z(a, b) if
 *   z < (a + 1) / (a + b + 2), and for I_{1-z}(b, a) if not.
 *
 * The terms of the power series shrink by a factor of about z per term
 * (1 - z for the series of I_{1-z}(b, a)), so the series is limited to
 * z <= 0.5. The condition b z <= 1 does not bound z by itself when b < 1.
 * Without the limit, the series would take most of its 600 terms at z = 0.95,
 * where the continued fraction needs tens of steps.
 *
 * The continued fraction converges rapidly for z < (a + 1) / (a + b + 2)
 * (Numerical Recipes, 3rd ed., section 6.4). When z is at or above that
 * point, 1 - z is at or below (b + 1) / (a + b + 2), the matching point for
 * the fraction of I_{1-z}(b, a). The mode (a - 1) / (a + b - 2) is not used
 * as the threshold because it exceeds 1 for b < 1. The fraction of I_z(a, b)
 * would then be used for z near 1, where it needs hundreds of steps.
 *
 * The partials are 0 at z = 0 and z = 1.
 *
 * @tparam grads_requested <code>G1</code> for the derivative with respect to
 *   a, <code>G2</code> for the derivative with respect to b,
 *   <code>BOTH</code> for both
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of z
 * @tparam T4 type of digammaA
 * @tparam T5 type of digammaB
 * @tparam T6 type of digammaSum
 * @tparam T7 type of betaAB
 * @param a a
 * @param b b
 * @param z z
 * @param digammaA the value of <code>digamma(a)</code>, not used for
 *   <code>G2</code>
 * @param digammaB the value of <code>digamma(b)</code>, not used for
 *   <code>G1</code>
 * @param digammaSum the value of <code>digamma(a + b)</code>
 * @param betaAB accepted for interface compatibility and not used. The
 *   normalization is evaluated in log space.
 * @return d/da I_z(a, b) for <code>G1</code>, d/db for <code>G2</code>, and
 *   the pair (d/da, d/db) for <code>BOTH</code>
 */
template <GRAD_RIB grads_requested, typename T1, typename T2, typename T3,
          typename T4, typename T5, typename T6, typename T7>
inline auto grad_reg_inc_beta(const T1& a, const T2& b, const T3& z,
                              const T4& digammaA, const T5& digammaB,
                              const T6& digammaSum, const T7& /* betaAB */) {
  // T7 is left out so that every branch returns the same type
  using ret_t = return_type_t<T1, T2, T3, T4, T5, T6>;
  // the partials requested after the symmetry relation swaps a and b
  constexpr GRAD_RIB swapped_grads
      = grads_requested == GRAD_RIB::G1
            ? GRAD_RIB::G2
            : grads_requested == GRAD_RIB::G2 ? GRAD_RIB::G1 : GRAD_RIB::BOTH;

  const double z_d = value_of_rec(z);
  if (z_d == 0.0 || z_d == 1.0) {
    if constexpr (grads_requested == GRAD_RIB::BOTH) {
      return std::pair<ret_t, ret_t>{0.0, 0.0};
    } else {
      return ret_t{0.0};
    }
  }

  const double a_d = value_of_rec(a);
  const double b_d = value_of_rec(b);
  const double mean = a_d / (a_d + b_d);
  const bool use_power_series
      = (z_d < mean && b_d * z_d <= 1.0 && z_d <= 0.5)
        || (z_d >= mean && a_d * (1.0 - z_d) <= 1.0 && z_d >= 0.5);
  if (use_power_series) {
    if (z_d >= mean) {
      return internal::inc_beta_power_series_kernel<swapped_grads, true>(
          b, a, 1.0 - z, digammaB, digammaA, digammaSum);
    } else {
      return internal::inc_beta_power_series_kernel<grads_requested, false>(
          a, b, z, digammaA, digammaB, digammaSum);
    }
  } else {
    if (z_d >= (a_d + 1.0) / (a_d + b_d + 2.0)) {
      return internal::inc_beta_continued_fraction_kernel<swapped_grads, true>(
          b, a, 1.0 - z, digammaB, digammaA, digammaSum);
    } else {
      return internal::inc_beta_continued_fraction_kernel<grads_requested,
                                                          false>(
          a, b, z, digammaA, digammaB, digammaSum);
    }
  }
}

/**
 * Computes the gradients of the regularized incomplete beta
 * function.  Specifically, this function computes gradients of
 * <code>inc_beta(a, b, z)</code> with respect to the arguments
 * <code>a</code> and <code>b</code>.
 *
 * Writes both partials. To compute only one, use the overload templated on
 * <code>GRAD_RIB</code>.
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
 * @param[in] betaAB the value of <code>beta(a, b)</code>, accepted for
 * interface compatibility and not used. The normalization is evaluated in
 * log space.
 */
template <typename T>
inline void grad_reg_inc_beta(T& g1, T& g2, const T& a, const T& b, const T& z,
                              const T& digammaA, const T& digammaB,
                              const T& digammaSum, const T& betaAB) {
  auto grads = grad_reg_inc_beta<GRAD_RIB::BOTH>(a, b, z, digammaA, digammaB,
                                                 digammaSum, betaAB);
  g1 = grads.first;
  g2 = grads.second;
}

}  // namespace math
}  // namespace stan
#endif
