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
 * A modified sign function that returns 1 for 0 inputs
 *
 * @tparam T Input type, can be either Eigen or scalar types
 * @param x Input argument to determine sign for
 * @return -1 for negative arguments otherwise 1
 */
template <typename T>
auto sign_binary(const T& x) {
  return sign(as_array_or_scalar(2 * sign(value_of_rec(x))) + 1);
}

/**
 * A modified log_sum_exp function that tracks the signs of the
 * input arguments and returns the appropriate sign for the output
 *
 * @tparam T Input type, can be either Eigen or scalar types
 * @param x Input argument to determine sign for
 * @return -1 for negative arguments otherwise 1
 */
template <typename T1>
inline value_type_t<T1> sum_exp_signed(const T1& v) {
  const value_type_t<T1> max_val = v.maxCoeff();
  const value_type_t<T1> val = sum(exp(v.array().abs() - max_val) * sign_binary(v));
  return exp(max_val + log(fabs(val))) * sign_binary(val);
}

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
                   double outer_precision, double inner_precision,
                   int outer_steps, int inner_steps) {
  using std::max;
  using scalar_t = return_type_t<Ta, Tb, Tz>;
  using MapA = Eigen::Map<Eigen::Matrix<value_type_t<Ta>, -1, 1>>;
  using MapB = Eigen::Map<Eigen::Matrix<value_type_t<Tb>, -1, 1>>;
  using MapT
      = Eigen::Map<Eigen::Matrix<scalar_t, -1, 1>, 0, Eigen::InnerStride<>>;
  using MapSignT = Eigen::Map<Eigen::VectorXi, 0, Eigen::InnerStride<>>;
  using Ta_plain = plain_type_t<Ta>;
  using Tb_plain = plain_type_t<Tb>;
  using T_vec = Eigen::Matrix<scalar_t, -1, 1>;
  ref_type_t<Ta> a_ref = a;
  ref_type_t<Tb> b_ref = b;
  int a_size = a.size();
  int b_size = b.size();

  bool condition_1 = (a_ref.size() > (b_ref.size() + 1)) && (z != 0);
  bool condition_2 = (a_ref.size() == (b_ref.size() + 1)) && (fabs(z) > 1);

  if (condition_1 || condition_2) {
    std::stringstream msg;
    msg << "hypergeometric pFq gradient does not meet convergence "
        << "conditions with given arguments. "
        << "a: " << a_ref << ", b: " << b_ref << ", "
        << "z: " << z;
    throw std::domain_error(msg.str());
  }

  Eigen::VectorXi a_signs = sign_binary(a_ref);
  Eigen::VectorXi b_signs = sign_binary(b_ref);
  int z_sign = sign_binary(z);

  Ta_plain ap1 = (a.array() + 1).matrix();
  Tb_plain bp1 = (b.array() + 1).matrix();
  Ta_plain log_a = log(a_ref);
  Tb_plain log_b = log(b_ref);
  scalar_type_t<Tz> log_z = log(fabs(z));
  scalar_type_t<Ta> log_a_prod = sum(log_a);
  scalar_type_t<Tb> log_b_prod = sum(log_b);
  scalar_type_t<Ta> a_prod = prod(a_ref);
  scalar_type_t<Tb> b_prod = prod(b_ref);

  // Only need the infinite sum for partials wrt a & b
  if (calc_a || calc_b) {
    double log_outer_precision = log(outer_precision);
    double log_inner_precision = log(inner_precision);

    // Append results at each iteration to a vector so that they
    // can be accumulated once, rather than every iteration
    std::vector<scalar_t> da_infsum;
    std::vector<scalar_t> db_infsum;
    //da_infsum.reserve(500 * a_size);
    //db_infsum.reserve(500 * b_size);

    T_vec da_mn = T_vec::Constant(a.size(), NEGATIVE_INFTY);
    T_vec db_mn = T_vec::Constant(b.size(), NEGATIVE_INFTY);

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

    double log_phammer_1n;
    double log_phammer_2_mpn;
    double lgamma_np1;
    T_vec log_phammer_an(a.size());
    T_vec log_phammer_ap1_n(a.size());
    T_vec log_phammer_bp1_n(b.size());
    T_vec log_phammer_bn(b.size());
    T_vec log_phammer_ap1_mpn(a.size());
    T_vec log_phammer_bp1_mpn(b.size());

    Eigen::VectorXi log_phammer_an_sign(a.size());
    Eigen::VectorXi log_phammer_ap1n_sign(a.size());
    Eigen::VectorXi log_phammer_bp1n_sign(b.size());
    Eigen::VectorXi log_phammer_bn_sign(b.size());
    Eigen::VectorXi log_phammer_ap1mpn_sign(a.size());
    Eigen::VectorXi log_phammer_bp1mpn_sign(b.size());
    Eigen::VectorXi log_phammer_ap1m_sign = Eigen::VectorXi::Ones(a.size());
    Eigen::VectorXi log_phammer_bp1m_sign = Eigen::VectorXi::Ones(b.size());

    std::vector<int> n_steps;
    while ((m < 50) && (m < outer_steps)) {
      int n = 0;
      Tz log_z_mn = log_z_m;
      int z_pow_mn_sign = z_sign;
      scalar_t inner_diff = 0;
      lgamma_np1 = 0;

      Ta_plain ap1n = ap1;
      Tb_plain bp1n = bp1;
      Ta_plain an = a_ref;
      Tb_plain bn = b_ref;
      Ta_plain ap1mn = ap1m;
      Tb_plain bp1mn = bp1m;

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

      while ((inner_diff > log_inner_precision) && (n < inner_steps)) {
        // Numerator term
        scalar_t term1_mn = log_z_mn + sum(log_phammer_ap1_mpn)
                            + log_phammer_1m + log_phammer_1n;
        // Denominator term
        scalar_t term2_mn = lgamma_mp1 + lgamma_np1 + sum(log_phammer_bp1_mpn)
                            + log_phammer_2_mpn;
        int base_sign = z_pow_mn_sign * log_phammer_ap1mpn_sign.prod()
                        * log_phammer_bp1mpn_sign.prod();

        if (calc_a) {
          // Division (on log scale) for the a & b partials
          da_mn = (term1_mn + log_phammer_an.array())
                  - (term2_mn + log_phammer_ap1_n.array());
          // Determine signs of each element
          Eigen::VectorXi curr_signs_da
              = (base_sign * log_phammer_an_sign.array()
                 * log_phammer_ap1n_sign.array())
                    .matrix();

          // Append results and signs for iteration
          for (int i = 0; i < da_mn.size(); i++) {
            da_infsum.push_back(da_mn[i] * curr_signs_da[i]);
          }
        }

        if (calc_b) {
          db_mn = (term1_mn + log_phammer_bn.array())
                  - (term2_mn + log_phammer_bp1_n.array());
          Eigen::VectorXi curr_signs_db
              = (base_sign * log_phammer_bn_sign.array()
                 * log_phammer_bp1n_sign.array())
                    .matrix();

          for (int i = 0; i < db_mn.size(); i++) {
            db_infsum.push_back(db_mn[i] * curr_signs_db[i]);
          }
        }

        // Series convergence assessed by whether the maximum term is
        //   smaller than the specified criteria (precision)
        inner_diff = std::max(da_mn.maxCoeff(), db_mn.maxCoeff());

        // Increment the input arguments and rising factorials
        log_z_mn += log_z;
        Ta_plain ap1mn = (ap1n.array() + m).matrix();
        Tb_plain bp1mn = (bp1n.array() + m).matrix();
        log_phammer_1n += log1p(n);
        log_phammer_2_mpn += log(2 + m + n);
        log_phammer_ap1_n += log(stan::math::fabs(ap1n));
        log_phammer_bp1_n += log(stan::math::fabs(bp1n));
        log_phammer_an += log(stan::math::fabs(an));
        log_phammer_bn += log(stan::math::fabs(bn));
        log_phammer_ap1_mpn += log(stan::math::fabs(ap1mn));
        log_phammer_bp1_mpn += log(stan::math::fabs(bp1mn));

        z_pow_mn_sign *= z_sign;
        log_phammer_ap1n_sign.array() *= sign_binary(ap1n);
        log_phammer_bp1n_sign.array() *= sign_binary(bp1n);
        log_phammer_an_sign.array() *= sign_binary(an);
        log_phammer_bn_sign.array() *= sign_binary(bn);
        log_phammer_ap1mpn_sign.array() *= sign_binary(ap1mn);
        log_phammer_bp1mpn_sign.array() *= sign_binary(bp1mn);

        n += 1;
        lgamma_np1 += log(n);
        ap1n.array() += 1;
        bp1n.array() += 1;
        an.array() += 1;
        bn.array() += 1;
        ap1mn.array() += 1;
        bp1mn.array() += 1;
      }

      n_steps.push_back(n);
      n_iter = n;

      log_z_m += log_z;
      log_phammer_1m += log1p(m);
      log_phammer_2m += log(2 + m);
      log_phammer_ap1_m += log(stan::math::fabs(ap1m));
      log_phammer_ap1m_sign.array() *= sign_binary(ap1m).array();
      log_phammer_bp1_m += log(stan::math::fabs(bp1m));
      log_phammer_bp1m_sign.array() *= sign_binary(bp1m).array();

      m += 1;

      lgamma_mp1 += log(m);
      ap1m.array() += 1;
      bp1m.array() += 1;
    }

    if (m == outer_steps) {
      throw_domain_error("grad_pFq", "k (internal counter)", outer_steps,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    std::cout << m + sum(n_steps) << std::endl;

    std::cout << db_infsum.size() << std::endl;
    std::cout << da_infsum.size() << std::endl;

    if (calc_a) {
      T_vec da(a.size());
      for (int i = 0; i < a.size(); i++) {
        MapT iter_val(da_infsum.data() + i, da_infsum.size() / a.size(), Eigen::InnerStride<>(a.size()));
        da[i] = sum_exp_signed(iter_val);
      }

      auto pre_mult_a = (z * a_prod / a_ref.array() / b_prod).matrix();

      // Evaluate gradients into provided containers
      std::get<0>(grad_tuple) = pre_mult_a.cwiseProduct(da);
    }

    if (calc_b) {
      T_vec db(b.size());
      for (int i = 0; i < b.size(); i++) {
        MapT iter_val(db_infsum.data() + i, db_infsum.size(), Eigen::InnerStride<>(b.size()));
        db[i] = sum_exp_signed(iter_val);
      }
      auto pre_mult_b = ((z * a_prod) / (b.array() * b_prod)).matrix();
      std::get<1>(grad_tuple) = -pre_mult_b.cwiseProduct(db);
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple)
        = (a_prod / b_prod) * hypergeometric_pFq(ap1, bp1, z);
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
auto grad_pFq(const Ta& a, const Tb& b, const Tz& z,
              double outer_precision = 1e-9,
              double inner_precision = 1e-9,
              int outer_steps = 1e6,
              int inner_steps = 1e6) {
  using partials_t = partials_return_t<Ta, Tb, Tz>;
  std::tuple<promote_scalar_t<partials_t, plain_type_t<Ta>>,
             promote_scalar_t<partials_t, plain_type_t<Tb>>,
             promote_scalar_t<partials_t, plain_type_t<Tz>>>
      ret_tuple;
  internal::grad_pFq_impl<!is_constant<Ta>::value, !is_constant<Tb>::value,
                          !is_constant<Tz>::value>(
      ret_tuple, value_of(a), value_of(b), value_of(z), outer_precision,
      inner_precision, outer_steps, inner_steps);
  return ret_tuple;
}

}  // namespace math
}  // namespace stan
#endif
