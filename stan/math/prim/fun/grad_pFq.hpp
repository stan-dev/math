#ifndef STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP
#define STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_rising_factorial.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

/**
 * Returns the gradient of generalised hypergeometric function wrt to the
 * input arguments:
 * \f$ _pF_q(a_1,...,a_p;b_1,...,b_q;z) \f$
 *
 * The derivatives wrt a and b are defined using a Kampé de Fériet function,
 *   (https://en.wikipedia.org/wiki/Kamp%C3%A9_de_F%C3%A9riet_function).
 *   This is implemented below as an infinite sum (on the log scale) until
 *   convergence.
 *
 * \f$ \frac{\partial}{\partial a_1} =
 *     \frac{z\prod_{j=2}^pa_j}{\prod_{j=1}^qb_j}
 *     F_{q+1\:0\;1}^{p\:1\:2}\left(\begin{array}&a_1+1,...,a_p+1;1;1,a_1\\
 *        2, b_1+1,...,b_1+1;;a_1+1\end{array};z,z\right) \f$
 *
 * \f$ \frac{\partial}{\partial b_1}=
 *     -\frac{z\prod_{j=1}^pa_j}{b_1\prod_{j=1}^qb_j}
 *     F_{q+1\:0\;1}^{p\:1\:2}\left(\begin{array}&a_1+1,...,a_p+1;1;1,b_1\\
 *        2, b_1+1,...,b_1+1;;b_1+1\end{array};z,z\right) \f$
 *
 * \f$ \frac{\partial}{\partial z}= \frac{\prod_{j=1}^pa_j}{\prod_{j=1}^qb_j}
 *      {}_pF_q(a_1 + 1,...,a_p + 1; b_1 +1,...,b_q+1;z) \f$
 *
 * @tparam calc_a Boolean for whether to calculate derivatives wrt to 'a'
 * @tparam calc_b Boolean for whether to calculate derivatives wrt to 'b'
 * @tparam calc_z Boolean for whether to calculate derivatives wrt to 'z'
 * @tparam TupleT Type of tuple containing objects to evaluate gradients into
 * @tparam Ta Eigen type with either one row or column at compile time
 * @tparam Tb Eigen type with either one row or column at compile time
 * @tparam Tz Scalar type
 * @param[in] grad_tuple Tuple of references to evaluate gradients into
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @param[in] precision Convergence criteria for infinite sum
 * @param[in] max_steps Maximum number of iterations for infinite sum
 * @return Generalised hypergeometric function
 */
template <bool calc_a, bool calc_b, bool calc_z, typename TupleT, typename Ta,
          typename Tb, typename Tz,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
void grad_pFq_impl(TupleT&& grad_tuple, const Ta& a, const Tb& b, const Tz& z,
                   double precision, int max_steps) {
  using std::max;
  using scalar_t = return_type_t<Ta, Tb, Tz>;
  using MapT
      = Eigen::Map<Eigen::Matrix<scalar_t, -1, 1>, 0, Eigen::InnerStride<>>;
  using Ta_plain = plain_type_t<Ta>;
  using Tb_plain = plain_type_t<Tb>;
  using T_vec = Eigen::Matrix<scalar_t, -1, 1>;

  Ta_plain ap1 = (a.array() + 1).matrix();
  Tb_plain bp1 = (b.array() + 1).matrix();
  Ta_plain log_a = log(a);
  Tb_plain log_b = log(b);
  scalar_type_t<Tz> log_z = log(z);
  scalar_type_t<Ta> sum_log_a = sum(log_a);
  scalar_type_t<Tb> sum_log_b = sum(log_b);

  // Declare vectors to accumulate sums into
  // where NEGATIVE_INFTY is zero on the log scale
  T_vec da_infsum(a.size());
  std::vector<scalar_t> da_infsumt(0);
  T_vec db_infsum(b.size());
  std::vector<scalar_t> db_infsumt(0);
  // Vectors to accumulate outer sum into
  T_vec da_iter_m = T_vec::Constant(a.size(), NEGATIVE_INFTY);
  T_vec db_iter_m = T_vec::Constant(b.size(), NEGATIVE_INFTY);
  T_vec da_mn = T_vec::Constant(a.size(), NEGATIVE_INFTY);
  T_vec db_mn = T_vec::Constant(b.size(), NEGATIVE_INFTY);

  // Only need the infinite sum for partials wrt p & b
  if (calc_a || calc_b) {
    double outer_precision = log(precision);
    double inner_precision = log(precision) - LOG_TWO;

    int m = 0;
    scalar_t outer_diff = 0;

    double lgamma_mp1 = 0;
    double log_phammer_1m = 0;
    double log_phammer_2m = 0;
    T_vec log_phammer_ap1_m = T_vec::Zero(ap1.size());
    T_vec log_phammer_bp1_m = T_vec::Zero(bp1.size());

    double log_phammer_1n;
    double log_phammer_2_mpn;
    double lgamma_np1;
    T_vec log_phammer_an(a.size());
    T_vec log_phammer_ap1_n(a.size());
    T_vec log_phammer_bp1_n(b.size());
    T_vec log_phammer_bn(b.size());
    T_vec log_phammer_ap1_mpn(a.size());
    T_vec log_phammer_bp1_mpn(b.size());

    while ((outer_diff > outer_precision) && (m < max_steps)) {
      // Vectors to accumulate outer sum into
      std::vector<scalar_t> da_iter_mt(0);
      std::vector<scalar_t> db_iter_mt(0);

      int n = 0;
      scalar_t inner_diff = 0;
      lgamma_np1 = 0;

      log_phammer_1n = 0;
      log_phammer_an.setZero();
      log_phammer_ap1_n.setZero();
      log_phammer_bp1_n.setZero();
      log_phammer_bn.setZero();
      log_phammer_ap1_mpn = log_phammer_ap1_m;
      log_phammer_bp1_mpn = log_phammer_bp1_m;
      log_phammer_2_mpn = log_phammer_2m;

      while ((inner_diff > inner_precision) & (n < (max_steps / 2))) {
        // Numerator term
        scalar_t term1_mn = (m + n) * log_z + sum(log_phammer_ap1_mpn)
                            + log_phammer_1m + log_phammer_1n;
        // Denominator term
        scalar_t term2_mn = lgamma_mp1 + lgamma_np1 + sum(log_phammer_bp1_mpn)
                            + log_phammer_2_mpn;
        if (calc_a) {
          // Division (on log scale) for the a & b partials
          da_mn = (term1_mn + log_phammer_an.array())
                  - (term2_mn + log_phammer_ap1_n.array());

          // Perform a row-wise log_sum_exp to accumulate current iteration
          for (int i = 0; i < da_mn.size(); i++) {
            da_iter_mt.emplace_back(da_mn[i]);
          }
        }
        if (calc_b) {
          db_mn = (term1_mn + log_phammer_bn.array())
                  - (term2_mn + log_phammer_bp1_n.array());
          for (int i = 0; i < db_mn.size(); i++) {
            db_iter_mt.emplace_back(db_mn[i]);
          }
        }

        // Series convergence assessed by whether the sum of all terms is
        //   smaller than the specified criteria (precision)
        inner_diff = std::max(da_mn.maxCoeff(), db_mn.maxCoeff());

        log_phammer_1n += log1p(n);
        log_phammer_2_mpn += log1p(1 + m + n);
        log_phammer_ap1_n.array() += log(ap1.array() + n);
        log_phammer_bp1_n.array() += log(bp1.array() + n);
        log_phammer_an.array() += log(a.array() + n);
        log_phammer_bn.array() += log(b.array() + n);
        log_phammer_ap1_mpn.array() += log(ap1.array() + m + n);
        log_phammer_bp1_mpn.array() += log(bp1.array() + m + n);
        n += 1;
        lgamma_np1 += log(n);
      }
      if (calc_a) {
        for (int i = 0; i < a.size(); i++) {
          da_iter_m[i] = log_sum_exp(
              MapT(da_iter_mt.data() + i, n, Eigen::InnerStride<>(a.size())));
          da_infsumt.emplace_back(da_iter_m[i]);
        }
      }
      if (calc_b) {
        for (int i = 0; i < b.size(); i++) {
          db_iter_m[i] = log_sum_exp(
              MapT(db_iter_mt.data() + i, n, Eigen::InnerStride<>(b.size())));
          db_infsumt.emplace_back(db_iter_m[i]);
        }
      }

      // Assess convergence of outer loop
      outer_diff = std::max(da_iter_m.maxCoeff(), db_iter_m.maxCoeff());

      log_phammer_1m += log1p(m);
      log_phammer_2m += log1p(1 + m);
      log_phammer_ap1_m.array() += log(ap1.array() + m);
      log_phammer_bp1_m.array() += log(bp1.array() + m);

      m += 1;

      lgamma_mp1 += log(m);
    }
    if (m == max_steps) {
      throw_domain_error("grad_pFq", "k (internal counter)", max_steps,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    if (calc_a) {
      for (int i = 0; i < a.size(); i++) {
        da_infsum[i] = log_sum_exp(
            MapT(da_infsumt.data() + i, m, Eigen::InnerStride<>(a.size())));
      }
      // Workaround to construct vector where each element is the product of
      //   all other elements
      Eigen::VectorXi ind_vector
          = Eigen::VectorXi::LinSpaced(a.size(), 0, a.size() - 1);

      Ta_plain prod_excl_curr = ind_vector.unaryExpr(
          [&log_a, &sum_log_a](int i) { return sum_log_a - log_a[i]; });
      T_vec pre_mult_a = (log_z + prod_excl_curr.array() - sum_log_b).matrix();

      // Evaluate gradients into provided containers
      std::get<0>(grad_tuple) = exp(pre_mult_a + da_infsum);
    }

    if (calc_b) {
      for (int i = 0; i < b.size(); i++) {
        db_infsum[i] = log_sum_exp(
            MapT(db_infsumt.data() + i, m, Eigen::InnerStride<>(b.size())));
      }
      T_vec pre_mult_b = (log_z + sum_log_a) - (log_b.array() + sum_log_b);
      std::get<1>(grad_tuple) = -exp(pre_mult_b + db_infsum);
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple)
        = exp(sum_log_a - sum_log_b) * hypergeometric_pFq(ap1, bp1, z);
  }
}
}  // namespace internal

/**
 * Wrapper function for calculating gradients for the generalized
 * hypergeometric function. The function always returns a tuple with
 * three elements (gradients wrt a, b, and z, respectively), but the
 * elements will only be defined/calculated when the respective parameter
 * is not a primitive type.
 *
 * @tparam Ta Eigen type with either one row or column at compile time
 * @tparam Tb Eigen type with either one row or column at compile time
 * @tparam Tz Scalar type
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @param[in] precision Convergence criteria for infinite sum
 * @param[in] max_steps Maximum number of iterations for infinite sum
 * @return Tuple of gradients
 */
template <typename Ta, typename Tb, typename Tz>
auto grad_pFq(const Ta& a, const Tb& b, const Tz& z, double precision = 1e-10,
              int max_steps = 1e6) {
  using partials_t = partials_return_t<Ta, Tb, Tz>;
  std::tuple<promote_scalar_t<partials_t, plain_type_t<Ta>>,
             promote_scalar_t<partials_t, plain_type_t<Tb>>,
             promote_scalar_t<partials_t, plain_type_t<Tz>>>
      ret_tuple;
  internal::grad_pFq_impl<!is_constant<Ta>::value, !is_constant<Tb>::value,
                          !is_constant<Tz>::value>(
      ret_tuple, value_of(a), value_of(b), value_of(z), precision, max_steps);
  return ret_tuple;
}

}  // namespace math
}  // namespace stan
#endif
