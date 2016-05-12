#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BUILD_PRIM_BINARY_MATRIX_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BUILD_PRIM_BINARY_MATRIX_HPP

#include <vector>
#include <Eigen/Dense>

template <typename T, int R, int C>
static inline Eigen::Matrix<T, R, C>
build_prim_binary_matrix(T val, const Eigen::Matrix<T, R, C>& x) {
  using Eigen::Matrix;
  using std::vector;

  Matrix<T, R, C> result_matrix(x.rows(), x.cols());
  for (int i = 0; i < x.size(); ++i)
    result_matrix(i) = val;
  return result_matrix;
}
#endif
