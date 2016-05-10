#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_BINARY_MATRIX_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_BUILD_BINARY_MATRIX_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F, typename T, int R, int C>
static inline Eigen::Matrix<T, R, C> build_fwd_binary_matrix1(
const Eigen::Matrix<T, R, C>& x, int seed_index = -1) {
  using std::vector;

  Eigen::Matrix<T, R, C> result_matrix(x.rows(), x.cols());
  //Fills matrix with copies of valid_input values
  vector<T> inputs = build_binary_vector1<F>(vector<T>());
  for (int i = 0; i < x.size(); ++i) {
    if (seed_index == i)
      result_matrix(i) = build_binary_vector1<F>(vector<T>(), 
      seed_index % inputs.size())[i % inputs.size()];
    else
      result_matrix(i) = inputs[i % inputs.size()];
  }
  return result_matrix;
}

template <typename F, typename T, int R, int C>
static inline Eigen::Matrix<T, R, C> build_fwd_binary_matrix1(
int val_index, const Eigen::Matrix<T, R, C>& x, int seed_index = -1) {
  using std::vector;

  Eigen::Matrix<T, R, C> result_matrix(x.rows(), x.cols());
  vector<T> val_v = build_binary_vector1<F>(vector<T>());
  for (int i = 0; i < x.size(); ++i)
  {
    if (seed_index == i)
      result_matrix(i) = build_binary_vector1<F>(vector<T>(), 
      seed_index % val_v.size())[val_index];
    else
      result_matrix(i) = val_v[val_index];
  }
  return result_matrix;
}

template <typename F, typename T, int R, int C>
static inline Eigen::Matrix<T, R, C> build_fwd_binary_matrix2(
const Eigen::Matrix<T, R, C>& x, int seed_index = -1) {
  using std::vector;

  Eigen::Matrix<T, R, C> result_matrix(x.rows(), x.cols());
  //Fills matrix with copies of valid_input values
  vector<T> inputs = build_binary_vector2<F>(std::vector<T>());
  for (int i = 0; i < x.size(); ++i) {
    if (seed_index == i)
      result_matrix(i) = build_binary_vector2<F>(vector<T>(), 
      seed_index % inputs.size())[i % inputs.size()];
    else
      result_matrix(i) = inputs[i % inputs.size()];
  }
  return result_matrix;
}

template <typename F, typename T, int R, int C>
static inline Eigen::Matrix<T, R, C> build_fwd_binary_matrix2(
int val_index, const Eigen::Matrix<T, R, C>& x, int seed_index = -1) {
  using std::vector;

  Eigen::Matrix<T, R, C> result_matrix(x.rows(), x.cols());
  vector<T> val_v = build_binary_vector2<F>(vector<T>());
  for (int i = 0; i < x.size(); ++i)
  {
    if (seed_index == i)
      result_matrix(i) = build_binary_vector2<F>(vector<T>(), 
      seed_index % val_v.size())[val_index];
    else
      result_matrix(i) = val_v[val_index];
  }
  return result_matrix;
}

#endif
