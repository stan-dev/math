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
#include <stan/math/prim/fun/min.hpp>
#include <stan/math/prim/fun/log_sum_exp_signed.hpp>
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
template <bool calc_a, bool calc_b, bool calc_z, typename Ta,
          typename Tb, typename Tz,
          require_all_eigen_vector_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
auto grad_pFq_impl(const Ta& a, const Tb& b, const Tz& z,
                   double precision, int max_iter) {
  using std::max;
  using scalar_t = return_type_t<Ta, Tb, Tz>;
  using Ta_plain = plain_type_t<Ta>;
  using Tb_plain = plain_type_t<Tb>;
  using T_vec = Eigen::Matrix<scalar_t, -1, 1>;
  using ArrayTa = Eigen::Array<scalar_type_t<Ta>, -1, 1>;
  using ArrayTb = Eigen::Array<scalar_type_t<Tb>, -1, 1>;
  using ArrayInt = Eigen::Array<int, -1, 1>;

  ref_type_t<Ta> a_ref = a;
  ref_type_t<Tb> b_ref = b;
  std::tuple<promote_scalar_t<scalar_t, Ta_plain>,
              promote_scalar_t<scalar_t, Tb_plain>,
              scalar_t> grad_tuple;

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

  ArrayTa a_array = a_ref.array();
  ArrayTb b_array = b_ref.array();

  // Only need the infinite sum for partials wrt a & b
  if (calc_a || calc_b) {
    ArrayTa a_array_inf = a_array;
    ArrayTa log_phammer_a = log(a_array);
    ArrayTa digamma_a_term = inv(a_array);

    ArrayTb b_array_inf = b_array;
    ArrayTb log_phammer_b = log(b_array);
    ArrayTb digamma_b_term = inv(b_array);

    Tz log_z = log(z);
    Tz log_z_pow_k = log_z;

    double log_factorial_k = 0;

    ArrayTa da = ArrayTa::Constant(a.size(), NEGATIVE_INFTY);
    ArrayTb db = ArrayTb::Constant(b.size(), NEGATIVE_INFTY);
    ArrayTa da_k = ArrayTa::Zero(a.size());
    ArrayTb db_k = ArrayTb::Zero(b.size());

    ArrayInt da_signs = ArrayInt::Ones(a.size());
    ArrayInt db_signs = ArrayInt::Ones(b.size());

    int k = 1;
    int min_iter = 10;
    scalar_t inner_diff = 1;

    while (((inner_diff > log(precision)) && (k < max_iter))) {
      scalar_t base_sum = (sum(log_phammer_a) + log_z_pow_k)
                          - (sum(log_phammer_b) + log_factorial_k);

      std::cout << "k: " << k  << std::endl;
      std::cout << "  base: " << base_sum << std::endl;
      if (calc_a) {
        da_k = log(digamma_a_term) + base_sum;
        std::cout << "  da_k: " << da_k.matrix().transpose() << std::endl;
        for (int i = 0; i < a.size(); i++) {
          std::forward_as_tuple(da(i), da_signs(i))
            = log_sum_exp_signed(da(i), da_signs(i), da_k(i), 1);
        }
        inner_diff = math::min(1.0, da_k.maxCoeff());
      }

      if (calc_b) {
        db_k = log(digamma_b_term) + base_sum;
        std::cout << "  db_k: " << db_k.matrix().transpose() << std::endl;
        for (int i = 0; i < b.size(); i++) {
          std::forward_as_tuple(db(i), db_signs(i))
            = log_sum_exp_signed(db(i), db_signs(i), db_k(i), -1);
        }
        inner_diff = math::min(1.0, db_k.maxCoeff());
      }

      a_array_inf += 1;
      b_array_inf += 1;

      log_phammer_a += log(a_array_inf);
      log_phammer_b += log(b_array_inf);
      log_z_pow_k += log_z;
      digamma_a_term += inv(a_array_inf);
      digamma_b_term += inv(b_array_inf);

      k += 1;
      log_factorial_k += log(k);
    }

    if (k == max_iter) {
      throw_domain_error("grad_pFq", "k (internal counter)", max_iter,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    if (calc_a) {
      std::get<0>(grad_tuple) = (exp(da) * da_signs).matrix();
    }

    if (calc_b) {
      std::get<1>(grad_tuple) = (exp(db) * db_signs).matrix();
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple) = (prod(a_ref) / prod(b_ref))
      * hypergeometric_pFq((a_array + 1).matrix(), (b_array + 1).matrix(), z);
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
auto grad_pFq(const Ta& a, const Tb& b, const Tz& z, double precision = 1e-10,
              int max_iter = 1e6) {
  return internal::grad_pFq_impl<!is_constant<Ta>::value, !is_constant<Tb>::value,
                          !is_constant<Tz>::value>(
      value_of(a), value_of(b), value_of(z), precision, max_iter);
}

}  // namespace math
}  // namespace stan
#endif
