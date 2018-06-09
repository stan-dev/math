#ifndef STAN_MATH_REV_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_REV_MAT_FUN_MATRIX_EXP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/prim/mat/fun/matrix_exp.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/rev/mat/fun/matrix_exp_action.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Overload matrix_exp function to perform matrix_exp_action.
 * This simplifies stan UI: we only expose
 * matrix_exp_action through the following signatures
 * - matrix_exp(A, B)
 * - matrix_exp(A, B, t)
 *
 * @tparam Ta scalar type matrix A
 * @tparam N Rows and cols matrix A, also rows of matrix B
 * @tparam Tb scalar type matrix B
 * @tparam Cb Columns matrix B
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplies B
 */
template <typename Ta, int N, typename Tb, int Cb>
inline typename boost::enable_if_c<boost::is_same<Ta, var>::value
                                       || boost::is_same<Tb, var>::value,
                                   Eigen::Matrix<var, N, Cb> >::type
inline matrix_exp(const Eigen::Matrix<Ta, N, N>& A,
                  const Eigen::Matrix<Tb, N, Cb>& B,
                  const double& t = 1.0) {
  return matrix_exp_action(A, B, t);
}

}  // namespace math
}  // namespace stan

#endif
