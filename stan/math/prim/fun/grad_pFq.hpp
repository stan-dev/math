#ifndef STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP
#define STAN_MATH_PRIM_FUN_GRAD_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/floor.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Returns the gradient of generalized hypergeometric function wrt to the
 * input arguments:
 * \f$ _pF_q(a_1,...,a_p;b_1,...,b_q;z) \f$
 *
 * Where:
 * \f$
 * \frac{\partial }{\partial a_1} =
 *  \sum_{k=0}^{\infty}{
 *    \frac
 *      {\psi\left(k+a_1\right)\left(\prod_{j=1}^p\left(a_j\right)_k\right)z^k}
 *      {k!\prod_{j=1}^q\left(b_j\right)_k}}
 *    - \psi\left(a_1\right){}_pF_q(a_1,...,a_p;b_1,...,b_q;z)
 * \f$
 * \f$
 * \frac{\partial }{\partial b_1} =
 *  \psi\left(b_1\right){}_pF_q(a_1,...,a_p;b_1,...,b_q;z) -
 *  \sum_{k=0}^{\infty}{
 *    \frac
 *      {\psi\left(k+b_1\right)\left(\prod_{j=1}^p\left(a_j\right)_k\right)z^k}
 *      {k!\prod_{j=1}^q\left(b_j\right)_k}}
 * \f$
 *
 * \f$
 *  \frac{\partial }{\partial z} =
 *  \frac{\prod_{j=1}^{p}(a_j)}{\prod_{j=1}^{q} (b_j)}\
 *    {}_pF_q(a_1+1,...,a_p+1;b_1+1,...,b_q+1;z)
 * \f$
 *
 * Noting the the recurrence relation for the digamma function:
 * \f$ \psi(x + 1) = \psi(x) + \frac{1}{x} \f$, the gradients for the
 * function with respect to a and b then simplify to:
 * \f$
 * \frac{\partial }{\partial a_1} =
 *  \sum_{k=1}^{\infty}{
 *    \frac
 *      {\left(1 + \sum_{m=0}^{k-1}\frac{1}{m+a_1}\right)
 *        * \left(\prod_{j=1}^p\left(a_j\right)_k\right)z^k}
 *      {k!\prod_{j=1}^q\left(b_j\right)_k}}
 *    - {}_pF_q(a_1,...,a_p;b_1,...,b_q;z)
 * \f$
 * \f$
 * \frac{\partial }{\partial b_1} =
 *  {}_pF_q(a_1,...,a_p;b_1,...,b_q;z) -
 *  \sum_{k=1}^{\infty}{
 *    \frac
 *      {\left(1 + \sum_{m=0}^{k-1}\frac{1}{m+b_1}\right)
 *        * \left(\prod_{j=1}^p\left(a_j\right)_k\right)z^k}
 *      {k!\prod_{j=1}^q\left(b_j\right)_k}}
 * \f$
 *
 * @tparam calc_a Boolean for whether to calculate derivatives wrt to 'a'
 * @tparam calc_b Boolean for whether to calculate derivatives wrt to 'b'
 * @tparam calc_z Boolean for whether to calculate derivatives wrt to 'z'
 * @tparam TpFq Scalar type of hypergeometric_pFq return value
 * @tparam Ta Eigen type with either one row or column at compile time
 * @tparam Tb Eigen type with either one row or column at compile time
 * @tparam Tz Scalar type
 * @param[in] pfq_val Return value from the hypergeometric_pFq function
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @param[in] precision Convergence criteria for infinite sum
 * @param[in] max_steps Maximum number of iterations
 * @return Tuple of gradients
 */
template <bool calc_a = true, bool calc_b = true, bool calc_z = true,
          typename TpFq, typename Ta, typename Tb, typename Tz,
          typename T_Rtn = return_type_t<Ta, Tb, Tz>,
          typename Ta_Rtn = promote_scalar_t<T_Rtn, plain_type_t<Ta>>,
          typename Tb_Rtn = promote_scalar_t<T_Rtn, plain_type_t<Tb>>>
std::tuple<Ta_Rtn, Tb_Rtn, T_Rtn> grad_pFq(const TpFq& pfq_val, const Ta& a,
                                           const Tb& b, const Tz& z,
                                           double precision = 1e-14,
                                           int max_steps = 1e6) {
  using std::max;
  using Ta_Array = Eigen::Array<return_type_t<Ta>, -1, 1>;
  using Tb_Array = Eigen::Array<return_type_t<Tb>, -1, 1>;

  Ta_Array a_array = as_column_vector_or_scalar(a).array();
  Tb_Array b_array = as_column_vector_or_scalar(b).array();

  std::tuple<Ta_Rtn, Tb_Rtn, T_Rtn> ret_tuple;

  if (calc_a || calc_b) {
    std::get<0>(ret_tuple).setConstant(a.size(), -pfq_val);
    std::get<1>(ret_tuple).setConstant(b.size(), pfq_val);
    Eigen::Array<T_Rtn, -1, 1> a_grad(a.size());
    Eigen::Array<T_Rtn, -1, 1> b_grad(b.size());

    int k = 0;
    int base_sign = 1;

    static constexpr double dbl_min = std::numeric_limits<double>::min();
    Ta_Array a_k = (a_array == 0.0).select(dbl_min, abs(a_array));
    Tb_Array b_k = (b_array == 0.0).select(dbl_min, abs(b_array));
    Tz log_z = log(abs(z));

    // Identify the number of iterations to needed for each element to sign
    // flip from negative to positive - rather than checking at each iteration
    Eigen::ArrayXi a_pos_k = (a_array < 0.0)
                                 .select(-floor(value_of_rec(a_array)), 0)
                                 .template cast<int>();
    int all_a_pos_k = a_pos_k.maxCoeff();
    Eigen::ArrayXi b_pos_k = (b_array < 0.0)
                                 .select(-floor(value_of_rec(b_array)), 0)
                                 .template cast<int>();
    int all_b_pos_k = b_pos_k.maxCoeff();
    Eigen::ArrayXi a_sign = select(a_pos_k == 0, 1, -1);
    Eigen::ArrayXi b_sign = select(b_pos_k == 0, 1, -1);

    int z_sign = sign(value_of_rec(z));

    Ta_Array digamma_a = Ta_Array::Ones(a.size());
    Tb_Array digamma_b = Tb_Array::Ones(b.size());

    T_Rtn curr_log_prec = NEGATIVE_INFTY;
    T_Rtn log_base = 0;
    while ((k < 10 || curr_log_prec > log(precision)) && (k <= max_steps)) {
      curr_log_prec = NEGATIVE_INFTY;
      if (calc_a) {
        a_grad = log(abs(digamma_a)) + log_base;
        std::get<0>(ret_tuple).array()
            += exp(a_grad) * base_sign * sign(value_of_rec(digamma_a));

        curr_log_prec = max(curr_log_prec, a_grad.maxCoeff());
        digamma_a += inv(a_k) * a_sign;
      }

      if (calc_b) {
        b_grad = log(abs(digamma_b)) + log_base;
        std::get<1>(ret_tuple).array()
            -= exp(b_grad) * base_sign * sign(value_of_rec(digamma_b));

        curr_log_prec = max(curr_log_prec, b_grad.maxCoeff());
        digamma_b += inv(b_k) * b_sign;
      }

      log_base += (sum(log(a_k)) + log_z) - (sum(log(b_k)) + log1p(k));
      base_sign *= z_sign * a_sign.prod() * b_sign.prod();

      // Wrap negative value handling in a conditional on iteration number so
      // branch prediction likely to ignore once positive
      if (k < all_a_pos_k) {
        // Avoid log(0) and 1/0 in next iteration by using smallest double
        //  - This is smaller than EPSILON, so the following iteration will
        //    still be 1.0
        a_k = (a_k == 1.0 && a_sign == -1)
                  .select(dbl_min, (a_k < 1.0 && a_sign == -1)
                                       .select(1.0 - a_k, a_k + 1.0 * a_sign));
        a_sign = select(k == a_pos_k - 1, 1, a_sign);
      } else {
        a_k += 1.0;

        if (k == all_a_pos_k) {
          a_sign.setOnes();
        }
      }

      if (k < all_a_pos_k) {
        b_k = (b_k == 1.0 && b_sign == -1)
                  .select(dbl_min, (b_k < 1.0 && b_sign == -1)
                                       .select(1.0 - b_k, b_k + 1.0 * b_sign));
        b_sign = select(k == b_pos_k - 1, 1, b_sign);
      } else {
        b_k += 1.0;

        if (k == all_b_pos_k) {
          b_sign.setOnes();
        }
      }

      k += 1;
    }
  }
  if (calc_z) {
    std::get<2>(ret_tuple)
        = hypergeometric_pFq(add(a, 1.0), add(b, 1.0), z) * prod(a) / prod(b);
  }
  return ret_tuple;
}

}  // namespace math
}  // namespace stan
#endif
