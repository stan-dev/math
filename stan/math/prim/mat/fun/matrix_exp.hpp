#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_HPP

#include <stan/math/prim/mat/fun/matrix_exp_pade.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_2x2.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the matrix exponential of the input
 * matrix.
 *
 * @tparam T type of elements in the matrix
 * @param[in] A Matrix to exponentiate.
 * @return Matrix exponential, dynamically-sized.
 * @throw <code>std::invalid_argument</code> if the input matrix
 * is not square.
 */
template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_exp(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A) {
  check_square("matrix_exp", "input matrix", A);
  if (A.size() == 0)
    return {};

  return (A.cols() == 2
          && square(value_of(A(0, 0)) - value_of(A(1, 1)))
                     + 4 * value_of(A(0, 1)) * value_of(A(1, 0))
                 > 0)
             ? matrix_exp_2x2(A)
             : matrix_exp_pade(A);
}

/**
 * Return the matrix exponential of the input
 * statically-sized square matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam N size of the input square matrix.
 * @param[in] A Matrix to exponentiate.
 * @return Matrix exponential, statically-sized.
 */
template <typename T, int N>
inline Eigen::Matrix<T, N, N> matrix_exp(const Eigen::Matrix<T, N, N>& A) {
  if (N == 0)
    return {};

  return (N == 2
          && square(value_of(A(0, 0)) - value_of(A(1, 1)))
                     + 4 * value_of(A(0, 1)) * value_of(A(1, 0))
                 > 0)
             ? matrix_exp_2x2(A)
             : matrix_exp_pade(A);
}

/**
 * Return the exponential of the input scalar when it's in
 * the form of Eigen matrix.
 *
 * @tparam T type of elements in the matrix
 * @return 1x1 Matrix exponential, statically-sized.
 */
template <typename T>
inline Eigen::Matrix<T, 1, 1> matrix_exp(const Eigen::Matrix<T, 1, 1>& A) {
  using std::exp;
  Eigen::Matrix<T, 1, 1> res;
  res << exp(A(0));
  return res;
}

}  // namespace math
}  // namespace stan

#endif
