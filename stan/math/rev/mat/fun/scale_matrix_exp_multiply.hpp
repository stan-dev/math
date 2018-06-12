#ifndef STAN_MATH_REV_MAT_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/rev/mat/fun/matrix_exp_multiply.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Wrapper of matrix_exp_action function for a more literal name.
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
scale_matrix_exp_multiply(const double& t, const Eigen::Matrix<Ta, N, N>& A,
                          const Eigen::Matrix<Tb, N, Cb>& B) {
  return matrix_exp_action(A, B, t);
}

}  // namespace math
}  // namespace stan

#endif
