#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_FWD_MATRIX_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_FWD_MATRIX_HPP

#include <stan/math/fwd/mat.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_vector.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename T, int R, int C>
static inline Eigen::Matrix<T, R, C> build_fwd_matrix(
    const Eigen::Matrix<T, R, C>& x, int seed_index = -1) {
  using Eigen::Matrix;
  using std::vector;

  Eigen::Matrix<T, R, C> fvar_matrix(x.rows(), x.cols());
  size_t num_inputs = F::valid_inputs().size();
  // Fills matrix with copies of valid_input values
  for (int i = 0; i < x.size(); ++i) {
    std::vector<T> inputs;
    if (seed_index == i)
      inputs = build_fwd_vector<F>(std::vector<T>(), seed_index % num_inputs);
    else
      inputs = build_fwd_vector<F>(std::vector<T>());
    fvar_matrix(i) = inputs[i % num_inputs];
  }
  return fvar_matrix;
}

#endif
