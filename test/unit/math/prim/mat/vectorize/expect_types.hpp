#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_TYPES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_TYPES_HPP

#include <test/unit/math/prim/mat/vectorize/expect_match_return_t.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F>
void expect_int_types() {
  using stan::test::expect_match_return_t;
  using std::vector;

  expect_match_return_t<F, double, int>();
  expect_match_return_t<F, vector<double>, vector<int> >();
  expect_match_return_t<F, vector<vector<double> >, vector<vector<int> > >();
}

template <typename F, typename T>
void expect_types() {
  using Eigen::Matrix;
  using stan::test::expect_match_return_t;
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> vector_t;
  typedef Eigen::Matrix<T, 1, Eigen::Dynamic> row_vector_t;

  expect_match_return_t<F, T, T>();
  expect_match_return_t<F, vector<T>, vector<T> >();
  expect_match_return_t<F, vector<vector<T> >, vector<vector<T> > >();
  expect_match_return_t<F, matrix_t, matrix_t>();
  expect_match_return_t<F, vector<matrix_t>, vector<matrix_t> >();
  expect_match_return_t<F, vector_t, vector_t>();
  expect_match_return_t<F, vector<vector_t>, vector<vector_t> >();
  expect_match_return_t<F, row_vector_t, row_vector_t>();
  expect_match_return_t<F, vector<row_vector_t>, vector<row_vector_t> >();
}

#endif
