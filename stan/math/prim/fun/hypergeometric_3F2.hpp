#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_3F2_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_3F2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/append_row.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {
namespace internal {
template <typename Ta, typename Tb, typename Tz,
          typename T_return = return_type_t<Ta, Tb, Tz>,
          typename ArrayAT = Eigen::Array<scalar_type_t<Ta>, 3, 1>,
          typename ArrayBT = Eigen::Array<scalar_type_t<Ta>, 3, 1>,
          require_all_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
T_return hypergeometric_3F2_infsum(const Ta& a, const Tb& b, const Tz& z,
                                   double precision = 1e-6,
                                   int max_steps = 1e5) {
  ArrayAT a_array = as_array_or_scalar(a);
  ArrayBT b_array = append_row(as_array_or_scalar(b), 1.0);
  check_3F2_converges("hypergeometric_3F2", a_array[0], a_array[1], a_array[2],
                      b_array[0], b_array[1], z);

  T_return t_acc = 1.0;
  T_return log_t = 0.0;
  T_return log_z = log(fabs(z));
  Eigen::ArrayXi a_signs = sign(value_of_rec(a_array));
  Eigen::ArrayXi b_signs = sign(value_of_rec(b_array));
  plain_type_t<decltype(a_array)> apk = a_array;
  plain_type_t<decltype(b_array)> bpk = b_array;
  int z_sign = sign(value_of_rec(z));
  int t_sign = z_sign * a_signs.prod() * b_signs.prod();

  int k = 0;
  while (k <= max_steps && log_t >= log(precision)) {
    // Replace zero values with 1 prior to taking the log so that we accumulate
    // 0.0 rather than -inf
    const auto& abs_apk = math::fabs((apk == 0).select(1.0, apk));
    const auto& abs_bpk = math::fabs((bpk == 0).select(1.0, bpk));
    T_return p = sum(log(abs_apk)) - sum(log(abs_bpk));
    if (p == NEGATIVE_INFTY) {
      return t_acc;
    }

    log_t += p + log_z;
    t_acc += t_sign * exp(log_t);

    if (is_inf(t_acc)) {
      throw_domain_error("hypergeometric_3F2", "sum (output)", t_acc,
                         "overflow hypergeometric function did not converge.");
    }
    k++;
    apk.array() += 1.0;
    bpk.array() += 1.0;
    a_signs = sign(value_of_rec(apk));
    b_signs = sign(value_of_rec(bpk));
    t_sign = a_signs.prod() * b_signs.prod() * t_sign;
  }
  if (k == max_steps) {
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
auto hypergeometric_3F2(const Ta& a, const Tb& b, const Tz& z) {
  check_3F2_converges("hypergeometric_3F2", a[0], a[1], a[2], b[0], b[1], z);
  // Boost's pFq throws convergence errors in some cases, fallback to naive
  // infinite-sum approach (tests pass for these)
  if (z == 1.0 && (sum(b) - sum(a)) < 0.0) {
    return internal::hypergeometric_3F2_infsum(a, b, z);
  }
  return hypergeometric_pFq(to_vector(a), to_vector(b), z);
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
auto hypergeometric_3F2(const std::initializer_list<Ta>& a,
                        const std::initializer_list<Tb>& b, const Tz& z) {
  return hypergeometric_3F2(std::vector<Ta>(a), std::vector<Tb>(b), z);
}

}  // namespace math
}  // namespace stan
#endif
