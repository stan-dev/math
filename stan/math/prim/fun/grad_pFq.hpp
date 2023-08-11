#ifndef STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP
#define STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/log_rising_factorial.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <cmath>
#include <iostream>

namespace stan {
namespace math {
namespace internal {
/**
 * Returns the gradient of generalized hypergeometric function wrt to the
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
 * @param[in] outer_steps Maximum number of iterations for infinite sum
 * @param[in] inner_steps Maximum number of iterations for infinite sum
 * @return Generalized hypergeometric function
 */
template <bool calc_a, bool calc_b, bool calc_z, typename TupleT, typename Ta,
          typename Tb, typename Tz,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
void grad_pFq_impl(TupleT&& grad_tuple, const Ta& a, const Tb& b, const Tz& z,
                   double precision, int outer_steps, int inner_steps) {
  using std::max;
  using scalar_t = return_type_t<Ta, Tb, Tz>;
  using Ta_plain = plain_type_t<Ta>;
  using Tb_plain = plain_type_t<Tb>;
  using T_vec = Eigen::Matrix<scalar_t, -1, 1>;
  ref_type_t<Ta> a_ref_in = a;
  ref_type_t<Tb> b_ref_in = b;

  // Replace any zero inputs with the smallest number representable, so that
  // taking log and aggregating does not return -Inf
  Ta_plain a_ref
      = (a_ref_in.array() == 0).select(EPSILON, a_ref_in.array()).matrix();
  Tb_plain b_ref
      = (b_ref_in.array() == 0).select(EPSILON, b_ref_in.array()).matrix();
  int a_size = a.size();
  int b_size = b.size();

  // Convergence criteria for the gradients are the same as for the
  // Hypergeometric pFq itself
  bool condition_1 = (a_size > (b_size + 1)) && (z != 0);
  bool condition_2 = (a_size == (b_size + 1)) && (fabs(z) > 1);

  if (condition_1 || condition_2) {
    std::stringstream msg;
    msg << "hypergeometric pFq gradient does not meet convergence "
        << "conditions with given arguments. "
        << "a: " << a_ref << ", b: " << b_ref << ", "
        << "z: " << z;
    throw std::domain_error(msg.str());
  }

  // As the gradients will be aggregating on the log scale, we will track the
  // the values and the signs separately - to avoid taking the log of a
  // negative input
  Eigen::VectorXi a_signs = sign(value_of_rec(a_ref));
  Eigen::VectorXi b_signs = sign(value_of_rec(b_ref));
  int z_sign = sign(value_of_rec(z));

  Ta_plain ap1 = (a.array() + 1).matrix();
  Tb_plain bp1 = (b.array() + 1).matrix();
  scalar_type_t<Tz> log_z = log(fabs(z));
  scalar_type_t<Ta> a_prod = prod(a_ref);
  scalar_type_t<Tb> b_prod = prod(b_ref);

  // Only need the infinite sum for partials wrt a & b
  if (calc_a || calc_b) {
    double log_precision = log(precision);

    T_vec da_mn = T_vec::Constant(a_size, NEGATIVE_INFTY);
    T_vec db_mn = T_vec::Constant(b_size, NEGATIVE_INFTY);
    T_vec da = T_vec::Constant(a_size, 0.0);
    T_vec db = T_vec::Constant(b_size, 0.0);

    int m = 0;
    int n_iter = 2;

    double lgamma_mp1 = 0;
    double log_phammer_1m = 0;
    double log_phammer_2m = 0;
    T_vec log_phammer_ap1_m = T_vec::Zero(ap1.size());
    T_vec log_phammer_bp1_m = T_vec::Zero(bp1.size());
    Tz log_z_m = 0;
    Ta_plain ap1m = ap1;
    Tb_plain bp1m = bp1;

    Ta_plain ap1n = ap1;
    Tb_plain bp1n = bp1;
    Ta_plain an = a_ref;
    Tb_plain bn = b_ref;
    Ta_plain ap1mn = ap1m;
    Tb_plain bp1mn = bp1m;

    double log_phammer_1n;
    double log_phammer_2_mpn;
    double lgamma_np1;
    T_vec log_phammer_an(a_size);
    T_vec log_phammer_ap1_n(a_size);
    T_vec log_phammer_bp1_n(b_size);
    T_vec log_phammer_bn(b_size);
    T_vec log_phammer_ap1_mpn(a_size);
    T_vec log_phammer_bp1_mpn(b_size);

    int z_pow_m_sign = 1;
    Eigen::VectorXi curr_signs_da(a_size);
    Eigen::VectorXi curr_signs_db(b_size);
    Eigen::VectorXi log_phammer_an_sign(a_size);
    Eigen::VectorXi log_phammer_ap1n_sign(a_size);
    Eigen::VectorXi log_phammer_bp1n_sign(b_size);
    Eigen::VectorXi log_phammer_bn_sign(b_size);
    Eigen::VectorXi log_phammer_ap1mpn_sign(a_size);
    Eigen::VectorXi log_phammer_bp1mpn_sign(b_size);
    Eigen::VectorXi log_phammer_ap1m_sign = Eigen::VectorXi::Ones(a_size);
    Eigen::VectorXi log_phammer_bp1m_sign = Eigen::VectorXi::Ones(b_size);

    // If the inner loop converges in 1 iteration, then the sum has coverged
    // and another iteration of the outer loop is not needed
    while ((n_iter > 1) && (m < outer_steps)) {
      ap1n = ap1;
      bp1n = bp1;
      an = a_ref;
      bn = b_ref;
      ap1mn = ap1m;
      bp1mn = bp1m;

      int n = 0;
      Tz log_z_mn = log_z_m;
      int z_pow_mn_sign = z_pow_m_sign;
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
      log_phammer_an_sign.setOnes();
      log_phammer_ap1n_sign.setOnes();
      log_phammer_bp1n_sign.setOnes();
      log_phammer_bn_sign.setOnes();
      log_phammer_ap1mpn_sign = log_phammer_ap1m_sign;
      log_phammer_bp1mpn_sign = log_phammer_bp1m_sign;

      while ((inner_diff > log_precision) && (n < inner_steps)) {
        // Numerator term
        scalar_t term1_mn = log_z_mn + sum(log_phammer_ap1_mpn) + log_phammer_1m
                            + log_phammer_1n;
        // Denominator term
        scalar_t term2_mn = lgamma_mp1 + lgamma_np1 + sum(log_phammer_bp1_mpn)
                            + log_phammer_2_mpn;
        int base_sign = z_pow_mn_sign * log_phammer_ap1mpn_sign.prod()
                        * log_phammer_bp1mpn_sign.prod();

        if (calc_a) {
          // Division (on log scale) for the a & b partials
          // Determine signs of each element
          curr_signs_da = (base_sign * log_phammer_an_sign.array()
                           * log_phammer_ap1n_sign.array())
                              .matrix();
          da_mn = (term1_mn + log_phammer_an.array())
                  - (term2_mn + log_phammer_ap1_n.array());

          // Aggregate the sums on the natural scale, so that the sign can be
          // applied before aggregation
          da += exp(da_mn).cwiseProduct(curr_signs_da);
        }

        if (calc_b) {
          curr_signs_db = (base_sign * log_phammer_bn_sign.array()
                           * log_phammer_bp1n_sign.array())
                              .matrix();
          db_mn = (term1_mn + log_phammer_bn.array())
                  - (term2_mn + log_phammer_bp1_n.array());
          db += exp(db_mn).cwiseProduct(curr_signs_db);
        }

        // Series convergence assessed by whether the maximum term is
        //   smaller than the specified criteria (precision)
        inner_diff = max(da_mn.maxCoeff(), db_mn.maxCoeff());

        // Increment the input arguments and rising factorials
        log_z_mn += log_z;
        log_phammer_1n += log1p(n);
        log_phammer_2_mpn += log(2 + m + n);

        log_phammer_ap1_n.array()
            += log(math::fabs((ap1n.array() == 0).select(1.0, ap1n.array())));
        log_phammer_bp1_n.array()
            += log(math::fabs((bp1n.array() == 0).select(1.0, bp1n.array())));
        log_phammer_an.array()
            += log(math::fabs((an.array() == 0).select(1.0, an.array())));
        log_phammer_bn.array()
            += log(math::fabs((bn.array() == 0).select(1.0, bn.array())));
        log_phammer_ap1_mpn.array()
            += log(math::fabs((ap1mn.array() == 0).select(1.0, ap1mn.array())));
        log_phammer_bp1_mpn.array()
            += log(math::fabs((bp1mn.array() == 0).select(1.0, bp1mn.array())));

        z_pow_mn_sign *= z_sign;
        log_phammer_ap1n_sign.array() *= sign(value_of_rec(ap1n)).array();
        log_phammer_bp1n_sign.array() *= sign(value_of_rec(bp1n)).array();
        log_phammer_an_sign.array() *= sign(value_of_rec(an)).array();
        log_phammer_bn_sign.array() *= sign(value_of_rec(bn)).array();
        log_phammer_ap1mpn_sign.array() *= sign(value_of_rec(ap1mn)).array();
        log_phammer_bp1mpn_sign.array() *= sign(value_of_rec(bp1mn)).array();

        n += 1;
        lgamma_np1 += log(n);
        ap1n.array() += 1;
        bp1n.array() += 1;
        an.array() += 1;
        bn.array() += 1;
        ap1mn.array() += 1;
        bp1mn.array() += 1;
      }

      z_pow_m_sign *= z_sign;

      n_iter = n;

      log_z_m += log_z;
      log_phammer_1m += log1p(m);
      log_phammer_2m += log(2 + m);
      log_phammer_ap1_m += log(math::fabs(ap1m));
      log_phammer_ap1m_sign.array() *= sign(value_of_rec(ap1m)).array();
      log_phammer_bp1_m += log(math::fabs(bp1m));
      log_phammer_bp1m_sign.array() *= sign(value_of_rec(bp1m)).array();

      m += 1;

      lgamma_mp1 += log(m);
      ap1m.array() += 1;
      bp1m.array() += 1;
      ap1mn.array() += 1;
      bp1mn.array() += 1;
    }

    if (m == outer_steps) {
      throw_domain_error("grad_pFq", "k (internal counter)", outer_steps,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    if (calc_a) {
      auto pre_mult_a = (z * a_prod / a_ref.array() / b_prod).matrix();
      std::get<0>(grad_tuple) = std::move(pre_mult_a.cwiseProduct(da));
    }

    if (calc_b) {
      auto pre_mult_b = ((z * a_prod) / (b.array() * b_prod)).matrix();
      std::get<1>(grad_tuple) = std::move(-pre_mult_b.cwiseProduct(db));
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple)
        = std::move((a_prod / b_prod) * hypergeometric_pFq(ap1, bp1, z));
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
 * @param[in] outer_steps Maximum number of iterations for outer_steps
 * @param[in] inner_steps Maximum number of iterations for inner_steps
 * @return Tuple of gradients
 */
template <typename Ta, typename Tb, typename Tz>
auto grad_pFq(const Ta& a, const Tb& b, const Tz& z, double precision = 1e-10,
              int outer_steps = 1e6, int inner_steps = 1e6) {
  using partials_t = partials_return_t<Ta, Tb, Tz>;
  std::tuple<promote_scalar_t<partials_t, plain_type_t<Ta>>,
             promote_scalar_t<partials_t, plain_type_t<Tb>>,
             promote_scalar_t<partials_t, plain_type_t<Tz>>>
      ret_tuple;
  internal::grad_pFq_impl<!is_constant<Ta>::value, !is_constant<Tb>::value,
                          !is_constant<Tz>::value>(
      ret_tuple, value_of(a), value_of(b), value_of(z), precision, outer_steps,
      inner_steps);
  return ret_tuple;
}

}  // namespace math
}  // namespace stan
#endif
