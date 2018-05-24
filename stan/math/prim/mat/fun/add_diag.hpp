#ifndef STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP
#define STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/fmin.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>

namespace stan {
namespace math {
/**
 * Returns a Matrix with a constant added along the main diagonal.
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a type of element to add along the diagonal
 *
 * @param mat a matrix of dimensions greater than 1, 1
 * @param to_add a constant to add along the diagonal
 * @return a matrix with a constant added along main diagonal
 */
template <typename T_m, typename T_a>
inline typename Eigen::Matrix<typename return_type<T_m, T_a>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
add_diag(const Eigen::Matrix<T_m, -1, -1> &mat, const T_a &to_add) {
  size_t length_diag;
  length_diag = fmin(mat.rows(), mat.cols());
  Eigen::Matrix<typename return_type<T_m, T_a>::type, -1, -1> out = mat;
  Eigen::Matrix<typename return_type<T_m, T_a>::type, 1, -1> temp_vec =
      mat.diagonal();

  for (size_t i = 0; i < length_diag; ++i)
    temp_vec[i] += to_add;

  out.diagonal() = temp_vec;
  return out;
}

/**
 * Returns a Matrix with a vector added along the main diagonal.
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a type of element to add along the diagonal
 *
 * @param mat a matrix of dimensions greater than 1, 1
 * @param to_add a vector of elements to be added along the main diagonal
 * @return a matrix with a constant added along main diagonal
 */
template <typename T_m, typename T_a>
inline typename Eigen::Matrix<typename return_type<T_m, T_a>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
add_diag(const Eigen::Matrix<T_m, -1, -1> &mat,
         const Eigen::Matrix<T_a, 1, -1> &to_add) {
  Eigen::Matrix<typename return_type<T_m, T_a>::type, -1, -1> out = mat;
  Eigen::Matrix<typename return_type<T_m, T_a>::type, 1, -1> temp_vec =
      mat.diagonal();

  out.diagonal() += to_add;
  return out;
}

/**
 * Returns a Matrix with a vector added along the main diagonal.
 *
 * @tparam T_m type of element in Eigen::Matrix
 * @tparam T_a type of element to add along the diagonal
 *
 * @param mat a matrix of dimensions greater than 1, 1
 * @param to_add a vector of elements to be added along the main diagonal
 * @return a matrix with a constant added along main diagonal
 */
template <typename T_m, typename T_a>
inline typename Eigen::Matrix<typename return_type<T_m, T_a>::type,
                              Eigen::Dynamic, Eigen::Dynamic>
add_diag(const Eigen::Matrix<T_m, -1, -1> &mat,
         const Eigen::Matrix<T_a, -1, 1> &to_add) {
  Eigen::Matrix<typename return_type<T_m, T_a>::type, -1, -1> out = mat;
  Eigen::Matrix<typename return_type<T_m, T_a>::type, 1, -1> temp_vec =
      mat.diagonal();

  out.diagonal() += to_add;
  return out;
}
}
}
#endif
