#ifndef STAN_MATH_REV_FUN_READ_COV_MATRIX_HPP
#define STAN_MATH_REV_FUN_READ_COV_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/read_corr_L.hpp>
#include <stan/math/rev/fun/read_cov_L.hpp>
#include <stan/math/rev/fun/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/rev/fun/rows_dot_product.hpp>
#include <stan/math/prim/fun/read_cov_matrix.hpp>

namespace stan {
namespace math {

/**
 * A generally worse alternative to call prior to evaluating the
 * density of an elliptical distribution
 *
 * @tparam T_CPCs type of CPCs vector (must be a `var_value<T>` where `T`
 *  inherits from EigenBase)
 * @tparam T_sds type of sds vector (must be a `var_value<T>` where `T`
 *  inherits from EigenBase)
 * @param CPCs on (-1, 1)
 * @param sds on (0, inf)
 * @param log_prob the log probability value to increment with the Jacobian
 * @return Covariance matrix for specified partial correlations.
 */
template <typename T_CPCs, typename T_sds,
          require_all_var_vector_t<T_CPCs, T_sds>* = nullptr>
var_value<Eigen::MatrixXd> read_cov_matrix(const T_CPCs& CPCs, const T_sds& sds,
                                           scalar_type_t<T_CPCs>& log_prob) {
  return multiply_lower_tri_self_transpose(read_cov_L(CPCs, sds, log_prob));
}

/**
 * Builds a covariance matrix from CPCs and standard deviations
 *
 * @tparam T_CPCs type of CPCs vector (must be a `var_value<T>` where `T`
 *  inherits from EigenBase)
 * @tparam T_sds type of sds vector (must be a `var_value<T>` where `T`
 *  inherits from EigenBase)
 * @param CPCs in (-1, 1)
 * @param sds in (0, inf)
 */
template <typename T_CPCs, typename T_sds,
          require_all_var_vector_t<T_CPCs, T_sds>* = nullptr>
inline var_value<Eigen::MatrixXd> read_cov_matrix(const T_CPCs& CPCs,
                                                  const T_sds& sds) {
  size_t K = sds.rows();

  var_value<Eigen::MatrixXd> corr_L = read_corr_L(CPCs, K);
  var_value<Eigen::MatrixXd> prod
      = sds.val().matrix().asDiagonal() * corr_L.val();

  reverse_pass_callback([sds, corr_L, prod]() mutable {
    corr_L.adj() += sds.val().matrix().asDiagonal() * prod.adj();
    sds.adj() += (prod.adj().cwiseProduct(corr_L.val())).rowwise().sum();
  });

  return multiply_lower_tri_self_transpose(prod);
}

}  // namespace math
}  // namespace stan

#endif
