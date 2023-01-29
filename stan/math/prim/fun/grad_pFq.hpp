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
template <bool calc_a, bool calc_b,
          typename Ta, typename Tb, typename Tz,
          typename ScalarT = return_type_t<Ta, Tb, Tz>,
          typename TupleT
            = std::tuple<promote_scalar_t<ScalarT, plain_type_t<Ta>>,
                         promote_scalar_t<ScalarT, plain_type_t<Tb>>,
                         promote_scalar_t<ScalarT, plain_type_t<Tz>>>,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
TupleT grad_pFq_ab_infsum(const Ta& a, const Tb& b, const Tz& z,
                   double precision, int max_steps) {
  TupleT grad_tuple;
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
        << "conditions with given arguments: \n"
        << "a: [" << to_row_vector(a) << "]\n"
        << "b: [" << to_row_vector(b) << "]\n"
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
  return grad_tuple;
}

template <bool calc_a, bool calc_b, bool calc_z, typename Ta,
          typename Tb, typename Tz,
          typename ScalarT = return_type_t<Ta, Tb, Tz>,
          typename TupleT
            = std::tuple<promote_scalar_t<ScalarT, plain_type_t<Ta>>,
                         promote_scalar_t<ScalarT, plain_type_t<Tb>>,
                         promote_scalar_t<ScalarT, plain_type_t<Tz>>>,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
TupleT grad_pFq_impl(const Ta& a, const Tb& b, const Tz& z,
                   double precision, int max_steps) {
  TupleT grad_tuple;
  if (calc_a || calc_b) {
    grad_tuple = grad_pFq_ab_infsum<calc_a, calc_b>(a, b, z, precision, max_steps);
  }

  if (calc_z) {
    std::get<2>(grad_tuple)
        = std::move((a.prod() / b.prod())
          * hypergeometric_pFq((a.array() + 1).matrix(), (b.array() + 1).matrix(), z));
  }

  return grad_tuple;
}

template <bool calc_a, bool calc_b, bool calc_z, typename Ta,
          typename Tb, typename Tz,
          typename ScalarT = return_type_t<Ta, Tb, Tz>,
          typename TupleT
            = std::tuple<promote_scalar_t<ScalarT, plain_type_t<Ta>>,
                         promote_scalar_t<ScalarT, plain_type_t<Tb>>,
                         promote_scalar_t<ScalarT, plain_type_t<Tz>>>,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
TupleT grad_2F1_impl(const Ta& a, const Tb& b, const Tz& z,
                   double precision, int max_steps) {
  TupleT grad_tuple;
  bool euler_transform = false;
  try {
    check_2F1_converges("hypergeometric_2F1", a[0], a[1], b[0], z);
  } catch (const std::exception& e) {
    // Apply Euler's hypergeometric transformation if function
    // will not converge with current arguments
    check_2F1_converges("hypergeometric_2F1 (euler transform)", b[0] - a[0], a[1], b[0],
                        z / (z - 1));
    euler_transform = true;
  }
  euler_transform = false;
  if (euler_transform) {
    promote_scalar_t<ScalarT, plain_type_t<Ta>> a_euler(2);
    a_euler << a[1], b[0] - a[0];
    ScalarT z_euler = z / (z - 1);
    if (calc_a || calc_b) {
      // 'b' gradients under Euler transform require gradients from 'a'
      constexpr bool calc_a_euler = calc_a || calc_b;
      TupleT grad_tuple_ab = grad_pFq_ab_infsum<calc_a_euler, calc_b>(
          a_euler, b, z_euler, precision, max_steps);

      auto pre_mult_ab = inv(pow(1.0 - z, a[1]));
      if (calc_a) {
        std::get<0>(grad_tuple)[0] = -pre_mult_ab * std::get<0>(grad_tuple_ab)[1];

        auto hyper_da2 = hypergeometric_2F1(a_euler[0], a[1], b[0], z_euler);
        std::get<0>(grad_tuple)[1]
            = -pre_mult_ab * hyper_da2 * log1m(z)
              + pre_mult_ab * std::get<0>(grad_tuple_ab)[0];
      }
      if (calc_b) {
        std::get<1>(grad_tuple)[0]
            = pre_mult_ab
              * (std::get<0>(grad_tuple_ab)[1] + std::get<1>(grad_tuple_ab)[0]);
      }
    }
    if (calc_z) {
      auto hyper1 = hypergeometric_2F1(a_euler[0], a_euler[1], b[0], z_euler);
      auto hyper2 = hypergeometric_2F1(1 + a[1], 1 - a[0] + b[0], 1 + b[0], z_euler);
      auto pre_mult = a[1] * pow(1 - z, -1 - a[1]);
      std::get<2>(grad_tuple)
          = a[1] * pow(1 - z, -1 - a[1]) * hyper1
            + (a[1] * (b[0] - a[0]) * pow(1 - z, -a[1])
               * (inv(z - 1) - z / square(z - 1)) * hyper2)
                  / b[0];
    }
  } else {
    grad_tuple = grad_pFq_impl<calc_a, calc_b, calc_z>(
      a, b, z, precision, max_steps
    );
  }

  return grad_tuple;
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
  if (a.size() == 2 && b.size() == 1) {
  return internal::grad_2F1_impl<!is_constant<Ta>::value, !is_constant<Tb>::value,
                          !is_constant<Tz>::value>(
      value_of(a), value_of(b), value_of(z), precision, max_steps);
  }
  return internal::grad_pFq_impl<!is_constant<Ta>::value, !is_constant<Tb>::value,
                          !is_constant<Tz>::value>(
      value_of(a), value_of(b), value_of(z), precision, max_steps);
}

}  // namespace math
}  // namespace stan
#endif
