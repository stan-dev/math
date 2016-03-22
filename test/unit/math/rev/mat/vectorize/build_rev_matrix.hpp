#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_REV_MATRIX_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_BUILD_REV_MATRIX_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, int R, int C>
static inline Eigen::Matrix<stan::math::var, R, C>
build_rev_matrix(const Eigen::Matrix<stan::math::var, R, C>& x) {
  using Eigen::Matrix;
  using std::vector;
  using stan::math::var;

  Matrix<var, R, C> var_matrix(x.rows(), x.cols());
  vector<double> inputs = F::valid_inputs();
  //Fills matrix with copies of valid_input values
  for (int i = 0; i < x.size(); ++i)
    var_matrix(i) = inputs[i % inputs.size()];
  return var_matrix;
}

#endif
