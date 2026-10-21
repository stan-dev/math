#ifndef STAN_MATH_PRIM_FUN_INC_BETA_DDA_HPP
#define STAN_MATH_PRIM_FUN_INC_BETA_DDA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>

namespace stan {
namespace math {

/**
 * Returns the partial derivative of the regularized
 * incomplete beta function, I_{z}(a, b) with respect to a.
 *
 * Evaluated by <code>grad_reg_inc_beta</code>, computing only the derivative
 * with respect to a.
 *
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of z
 * @tparam T4 type of digamma_a
 * @tparam T5 type of digamma_ab
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @param digamma_a value of digamma(a)
 * @param digamma_ab value of digamma(a + b)
 * @return partial derivative of the incomplete beta with respect to a
 *
 * @pre a >= 0
 * @pre b >= 0
 * @pre 0 <= z <= 1
 */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          require_all_stan_scalar_t<T1, T2, T3, T4, T5>* = nullptr>
inline auto inc_beta_dda(const T1& a, const T2& b, const T3& z,
                         const T4& digamma_a, const T5& digamma_ab) {
  return grad_reg_inc_beta<GRAD_RIB::G1>(a, b, z, digamma_a, 0.0, digamma_ab,
                                         0.0);
}

/**
 * Enables the vectorized application of the inc_beta_dda function, when
 *  any arguments are vectors. Scalar arguments are broadcast against the
 *  vectors, which must all have the same size.
 *
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of z
 * @tparam T4 type of digamma_a
 * @tparam T5 type of digamma_ab
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @param digamma_a value of digamma(a)
 * @param digamma_ab value of digamma(a + b)
 * @return column array of the partial derivatives of the incomplete beta
 *  with respect to a
 * @throw std::invalid_argument if the vector arguments differ in size
 */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          require_any_vector_t<T1, T2, T3, T4, T5>* = nullptr,
          require_all_not_var_matrix_t<T1, T2, T3, T4, T5>* = nullptr>
inline auto inc_beta_dda(const T1& a, const T2& b, const T3& z,
                         const T4& digamma_a, const T5& digamma_ab) {
  check_consistent_sizes("inc_beta_dda", "a", a, "b", b, "z", z, "digamma_a",
                         digamma_a, "digamma_ab", digamma_ab);
  using ret_t = return_type_t<T1, T2, T3, T4, T5>;
  scalar_seq_view<T1> a_vec(a);
  scalar_seq_view<T2> b_vec(b);
  scalar_seq_view<T3> z_vec(z);
  scalar_seq_view<T4> digamma_a_vec(digamma_a);
  scalar_seq_view<T5> digamma_ab_vec(digamma_ab);
  const size_t N = max_size(a, b, z, digamma_a, digamma_ab);
  Eigen::Array<ret_t, Eigen::Dynamic, 1> result(N);
  for (size_t i = 0; i < N; ++i) {
    result.coeffRef(i) = inc_beta_dda(a_vec[i], b_vec[i], z_vec[i],
                                      digamma_a_vec[i], digamma_ab_vec[i]);
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
