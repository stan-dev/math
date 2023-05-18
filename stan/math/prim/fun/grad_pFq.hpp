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
  ref_type_t<Ta> a_ref_in = a;
  ref_type_t<Tb> b_ref_in = b;
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
  int k = 0;

  // Only need the infinite sum for partials wrt a & b
  if (calc_a || calc_b) {
    while ((inner_diff > log_precision) && (k < max_iter)) {
      if (calc_a) {
      }

      if (calc_b) {
      }
    }

    if (k == max_iter) {
      throw_domain_error("grad_pFq", "k (internal counter)", max_iter,
                         "exceeded ",
                         " iterations, hypergeometric function gradient "
                         "did not converge.");
    }

    if (calc_a) {
      std::get<0>(grad_tuple) = ;
    }

    if (calc_b) {
      std::get<1>(grad_tuple) = ;
    }
  }

  if (calc_z) {
    std::get<2>(grad_tuple) =
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
