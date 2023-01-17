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
 * @param[in] outer_precision Convergence criteria for infinite sum
 * @param[in] inner_precision Convergence criteria for infinite sum
 * @param[in] outer_steps Maximum number of iterations for infinite sum
 * @param[in] inner_steps Maximum number of iterations for infinite sum
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
  using Ta_plain = plain_type_t<Ta>;
  using Tb_plain = plain_type_t<Tb>;
  using T_vec = Eigen::Matrix<scalar_t, -1, 1>;
  ref_type_t<Ta> a_ref = a;
  ref_type_t<Tb> b_ref = b;

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

  scalar_type_t<Tz> log_z = log(fabs(z));

  // Only need the infinite sum for partials wrt a & b
  if (calc_a || calc_b) {
    double log_precision = log(precision);
    scalar_t inner_diff = 0;
    int m = 0;

    T_vec da_m = T_vec::Constant(a_size, NEGATIVE_INFTY);
    T_vec db_m = T_vec::Constant(b_size, NEGATIVE_INFTY);

    T_vec da = T_vec::Constant(a_size, 0.0);
    T_vec db = T_vec::Constant(b_size, 0.0);
    Ta_plain apm = a_ref;
    Tb_plain bpm = b_ref;

    // Replace zero values with 1, so that they have no impact on the logged accumulation
    Ta_plain log_phammer_a_mp1 = log(math::fabs((apm.array() == 0).select(1.0, apm.array())));
    Tb_plain log_phammer_b_mp1 = log(math::fabs((bpm.array() == 0).select(1.0, bpm.array())));

    Ta_plain inv_apm_sum = inv(a_ref);
    Tb_plain inv_bpm_sum = inv(b_ref);
    Tz log_z_pow_mp1 = log_z;
    double log_factorial_mp1 = 0;

    Eigen::VectorXi log_phammer_a_mp1_sign = a_signs;
    Eigen::VectorXi log_phammer_b_mp1_sign = b_signs;
    int log_z_pow_mp1_sign = z_sign;
    Eigen::VectorXi curr_signs_da(a_size);
    Eigen::VectorXi curr_signs_db(b_size);

    while ((inner_diff > log_precision) && (m < max_steps)) {
      // Numerator term
      scalar_t base_term = (sum(log_phammer_a_mp1) - sum(log_phammer_b_mp1))
                            + (log_z_pow_mp1 - log_factorial_mp1);

      int base_sign = log_z_pow_mp1_sign * log_phammer_a_mp1_sign.prod()
                      * log_phammer_b_mp1_sign.prod();

      if (calc_a) {
        da_m = base_term + log(fabs(inv_apm_sum)).array();
        curr_signs_da = base_sign * sign(value_of_rec(inv_apm_sum));

        // Aggregate the sums on the natural scale, so that the sign can be
        // applied before aggregation
        da += exp(da_m).cwiseProduct(curr_signs_da);
      }

      if (calc_b) {
        db_m = base_term + log(fabs(inv_bpm_sum)).array();
        curr_signs_db = base_sign * sign(value_of_rec(inv_bpm_sum));
        db -= exp(db_m).cwiseProduct(curr_signs_db);
      }

      // Series convergence assessed by whether the maximum term is
      //   smaller than the specified criteria (precision)
      inner_diff = max(da_m.maxCoeff(), db_m.maxCoeff());
      m += 1;

      // Increment the input arguments and rising factorials
      log_z_pow_mp1 += log_z;
      apm.array() += 1;
      bpm.array() += 1;
      inv_apm_sum.array() += inv(apm).array().isInf().select(0.0, inv(apm).array());
      inv_bpm_sum.array() += inv(bpm).array().isInf().select(0.0, inv(bpm).array());
      log_factorial_mp1 += log1p(m);

      log_phammer_a_mp1.array()
          += log(math::fabs((apm.array() == 0).select(1.0, apm.array())));
      log_phammer_b_mp1.array()
          += log(math::fabs((bpm.array() == 0).select(1.0, bpm.array())));

      log_z_pow_mp1_sign *= z_sign;
      log_phammer_a_mp1_sign.array() *= sign(value_of_rec(apm)).array();
      log_phammer_b_mp1_sign.array() *= sign(value_of_rec(bpm)).array();
    }

    if (m == max_steps) {
      throw_domain_error("grad_pFq", "k (internal counter)", max_steps,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    if (calc_a) {
      std::get<0>(grad_tuple) = std::move(da);
    }

    if (calc_b) {
      std::get<1>(grad_tuple) = std::move(db);
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple)
        = std::move((a_ref.prod() / b_ref.prod())
          * hypergeometric_pFq((a_ref.array() + 1).matrix(), (b_ref.array() + 1).matrix(), z));
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
auto grad_pFq(const Ta& a, const Tb& b, const Tz& z, double precision = 1e-14,
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
