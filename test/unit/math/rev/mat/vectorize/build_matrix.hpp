#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_MATRIX_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_MATRIX_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename T>
static inline Eigen::Matrix<stan::math::var, T::RowsAtCompileTime,
                              T::ColsAtCompileTime>
build_matrix(const T& x) {
  Eigen::Matrix<stan::math::var, T::RowsAtCompileTime, 
    T::ColsAtCompileTime> var_matrix(x.rows(), x.cols());
  std::vector<double> inputs = F::valid_inputs();
  for (int i = 0; i < x.size(); ++i) {
      var_matrix(i) = inputs[(i % inputs.size())];
  }
  return var_matrix;
}

#endif
