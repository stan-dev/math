#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_3F2_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_3F2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/append_row.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/hypergeometric_3F2_tail_bound.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <array>
#include <cmath>

namespace stan {
namespace math {
namespace internal {
/**
 * Sum the power series of the hypergeometric function 3F2 term by term.
 *
 * `hypergeometric_3F2` calls this function at z = 1 with sum(b) <= sum(a).
 * There `check_3F2_converges` accepts only a terminating series (a
 * polynomial), so the sum stops at the first zero term. It stops earlier
 * when a bound of the sum of all remaining terms is below `precision` times
 * the absolute value of the partial sum.
 *
 * @tparam Ta type of Eigen/Std vector 'a' arguments
 * @tparam Tb type of Eigen/Std vector 'b' arguments
 * @tparam Tz type of z argument
 * @param[in] a numerator parameters
 * @param[in] b denominator parameters
 * @param[in] z argument
 * @param[in] precision relative precision of the sum. The default 1e-17 is
 *   below half an ulp of the sum.
 * @param[in] max_steps number of steps to take
 * @return the sum of the series
 * @throw std::domain_error if the sum overflows or needs more than
 *   max_steps steps
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_3F2_infsum(
    const Ta& a, const Tb& b, const Tz& z, double precision = 1e-17,
    int max_steps = 1e5) {
  using T_return = return_type_t<Ta, Tb, Tz>;
  Eigen::Array<scalar_type_t<Ta>, 3, 1> a_array = as_array_or_scalar(a);
  Eigen::Array<scalar_type_t<Tb>, 3, 1> b_array
      = append_row(as_array_or_scalar(b), 1.0);
  check_3F2_converges("hypergeometric_3F2", a_array[0], a_array[1], a_array[2],
                      b_array[0], b_array[1], z);

  T_return t_acc = 1.0;
  T_return log_t = 0.0;
  auto log_z = log(fabs(z));
  Eigen::Array<int, 3, 1> a_signs = sign(value_of_rec(a_array));
  Eigen::Array<int, 3, 1> b_signs = sign(value_of_rec(b_array));
  int z_sign = sign(value_of_rec(z));
  int t_sign = z_sign * a_signs.prod() * b_signs.prod();

  // For the bound of the remaining terms: the parameters, and the index of
  // the last term that can be nonzero (the first zero numerator ends the
  // series, otherwise the steps end the sum)
  const std::array<double, 3> a_val{value_of_rec(a_array[0]),
                                    value_of_rec(a_array[1]),
                                    value_of_rec(a_array[2])};
  const std::array<double, 2> b_val{value_of_rec(b_array[0]),
                                    value_of_rec(b_array[1])};
  const double abs_z = std::fabs(value_of_rec(z));
  double last = max_steps + 1.0;
  for (double a_i : a_val) {
    if (a_i <= 0.0 && a_i == std::floor(a_i)) {
      last = std::fmin(last, -a_i);
    }
  }

  int k = 0;
  double abs_term = 1.0;
  while (k <= max_steps) {
    // A numerator parameter that has reached zero makes this term and every
    // later term zero: the series is a polynomial and has ended. Without
    // this stop the sign below is 0 while the magnitude keeps growing, and
    // 0 * inf gives NaN.
    if ((value_of_rec(a_array) == 0.0).any()) {
      return t_acc;
    }
    // Stop when the last term t_k and a bound of the sum of the remaining
    // terms t_{k + 1}, ..., t_last are negligible against the partial sum.
    // The ratios of the remaining terms are r_j for j = k, ..., last - 1.
    const double abs_sum = std::fabs(value_of_rec(t_acc));
    if (abs_term <= precision * abs_sum) {
      const double ratio_bound
          = hypergeometric_3F2_ratio_bound(a_val, b_val, abs_z, k, last - 1);
      const double tail
          = abs_term
            * hypergeometric_3F2_tail_sums(ratio_bound, last - k).first;
      if (tail <= precision * abs_sum) {
        return t_acc;
      }
    }
    // Replace zero values with 1 prior to taking the log so that we accumulate
    // 0.0 rather than -inf
    const auto& abs_apk = math::fabs((a_array == 0).select(1.0, a_array));
    const auto& abs_bpk = math::fabs((b_array == 0).select(1.0, b_array));
    auto p = sum(log(abs_apk)) - sum(log(abs_bpk));
    if (p == NEGATIVE_INFTY) {
      return t_acc;
    }

    log_t += p + log_z;
    const auto term = exp(log_t);
    t_acc += t_sign * term;
    abs_term = value_of_rec(term);

    if (is_inf(t_acc)) {
      throw_domain_error("hypergeometric_3F2", "sum (output)", t_acc,
                         "overflow hypergeometric function did not converge.");
    }
    k++;
    a_array += 1.0;
    b_array += 1.0;
    a_signs = sign(value_of_rec(a_array));
    b_signs = sign(value_of_rec(b_array));
    t_sign = z_sign * a_signs.prod() * b_signs.prod() * t_sign;
  }
  // The loop ends with k = max_steps + 1 when the steps run out
  if (k > max_steps) {
    throw_domain_error("hypergeometric_3F2", "k (internal counter)", max_steps,
                       "exceeded  iterations, hypergeometric function did not ",
                       "converge.");
  }
  return t_acc;
}
}  // namespace internal

/**
 * Hypergeometric function (3F2).
 *
 * Function reference: http://dlmf.nist.gov/16.2
 *
 * \f[
 *   _3F_2 \left(
 *     \begin{matrix}a_1 a_2 a[2] \\ b_1 b_2\end{matrix}; z
 *     \right) = \sum_k=0^\infty
 * \frac{(a_1)_k(a_2)_k(a_3)_k}{(b_1)_k(b_2)_k}\frac{z^k}{k!} \f]
 *
 * Where $(a_1)_k$ is an upper shifted factorial.
 *
 * Calculate the hypergeometric function (3F2) as the power series
 * directly to within <code>precision</code> or until
 * <code>max_steps</code> terms.
 *
 * This function does not have a closed form but will converge if:
 *   - <code>|z|</code> is less than 1
 *   - <code>|z|</code> is equal to one and <code>b[0] + b[1] < a[0] + a[1] +
 * a[2]</code> This function is a rational polynomial if
 *   - <code>a[0]</code>, <code>a[1]</code>, or <code>a[2]</code> is a
 *     non-positive integer
 * This function can be treated as a rational polynomial if
 *   - <code>b[0]</code> or <code>b[1]</code> is a non-positive integer
 *     and the series is terminated prior to the final term.
 *
 * @tparam Ta type of Eigen/Std vector 'a' arguments
 * @tparam Tb type of Eigen/Std vector 'b' arguments
 * @tparam Tz type of z argument
 * @param[in] a Always called with a[1] > 1, a[2] <= 0
 * @param[in] b Always called with int b[0] < |a[2]|,  <= 1)
 * @param[in] z z (is always called with 1 from beta binomial cdfs)
 * @param[in] precision precision of the infinite sum. defaults to 1e-6
 * @param[in] max_steps number of steps to take. defaults to 1e5
 * @return The 3F2 generalized hypergeometric function applied to the
 *  arguments {a1, a2, a3}, {b1, b2}
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
inline auto hypergeometric_3F2(const Ta& a, const Tb& b, const Tz& z) {
  check_size_match("hypergeometric_3F2", "a", a.size(), "3", 3);
  check_size_match("hypergeometric_3F2", "b", b.size(), "2", 2);

  auto a_ref = to_vector(a);
  auto b_ref = to_vector(b);

  check_3F2_converges("hypergeometric_3F2", a_ref[0], a_ref[1], a_ref[2],
                      b_ref[0], b_ref[1], z);
  // Boost's pFq throws convergence errors in some cases, fallback to naive
  // infinite-sum approach (tests pass for these). At z = 1 Boost also throws
  // when sum(b) == sum(a), also for a terminating series.
  if (z == 1.0 && (sum(b_ref) - sum(a_ref)) <= 0.0) {
    return internal::hypergeometric_3F2_infsum(a_ref, b_ref, z);
  }
  return hypergeometric_pFq(a_ref, b_ref, z);
}

/**
 * Hypergeometric function (3F2).
 *
 * Overload for initializer_list inputs
 *
 * @tparam Ta type of scalar 'a' arguments
 * @tparam Tb type of scalar 'b' arguments
 * @tparam Tz type of z argument
 * @param[in] a Always called with a[1] > 1, a[2] <= 0
 * @param[in] b Always called with int b[0] < |a[2]|,  <= 1)
 * @param[in] z z (is always called with 1 from beta binomial cdfs)
 * @param[in] precision precision of the infinite sum. defaults to 1e-6
 * @param[in] max_steps number of steps to take. defaults to 1e5
 * @return The 3F2 generalized hypergeometric function applied to the
 *  arguments {a1, a2, a3}, {b1, b2}
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_stan_scalar_t<Ta, Tb, Tz>* = nullptr>
inline auto hypergeometric_3F2(const std::initializer_list<Ta>& a,
                               const std::initializer_list<Tb>& b,
                               const Tz& z) {
  return hypergeometric_3F2(std::vector<Ta>(a), std::vector<Tb>(b), z);
}

}  // namespace math
}  // namespace stan
#endif
