#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_TYPES_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_TYPES_HPP

#include <test/unit/math/prim/mat/vectorize/expect_binary_match_return_t.hpp>
#include <vector>
#include <Eigen/Dense>

template <typename F>
void expect_int_binary_types() {
  using std::vector;
  using stan::test::expect_binary_match_return_t;

  expect_binary_match_return_t<F, double, int, int>();
  expect_binary_match_return_t<F, vector<vector<double> >, 
                               vector<vector<int> >, 
                               vector<vector<int> > >();
}

template <typename F, typename T1, typename T2>
void expect_binary_types() {
  using stan::test::expect_binary_match_return_t;
  using std::vector;
  using Eigen::Matrix; 
  typedef Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> matrix_t1;
  typedef Eigen::Matrix<T1, Eigen::Dynamic, 1> vector_t1;
  typedef Eigen::Matrix<T1, 1, Eigen::Dynamic> row_vector_t1;
  typedef Eigen::Matrix<T2, Eigen::Dynamic, Eigen::Dynamic> matrix_t2;
  typedef Eigen::Matrix<T2, Eigen::Dynamic, 1> vector_t2;
  typedef Eigen::Matrix<T2, 1, Eigen::Dynamic> row_vector_t2;
  typedef typename boost::math::tools::promote_args<T1, T2>::type result_t;
  typedef Eigen::Matrix<result_t, Eigen::Dynamic, Eigen::Dynamic>
          result_matrix_t;
  typedef Eigen::Matrix<result_t, Eigen::Dynamic, 1> result_vector_t;
  typedef Eigen::Matrix<result_t, 1, Eigen::Dynamic> result_row_vector_t;

  expect_binary_match_return_t<F, result_t, T1, T2>();
  expect_binary_match_return_t<F, vector<result_t>, vector<T1>, 
                               vector<T2> >();
  expect_binary_match_return_t<F, vector<vector<result_t> >, 
                               vector<vector<T1> >, 
                               vector< vector<T2> > >();
  expect_binary_match_return_t<F, result_matrix_t, matrix_t1, matrix_t2>(); 
  expect_binary_match_return_t<F, vector<result_matrix_t>, 
                               vector<matrix_t1>, vector<matrix_t2> >(); 
  expect_binary_match_return_t<F, result_vector_t, vector_t1, vector_t2>(); 
  expect_binary_match_return_t<F, vector<result_vector_t>, 
                               vector<vector_t1>, vector<vector_t2> >(); 
  expect_binary_match_return_t<F, result_row_vector_t, row_vector_t1, 
                               row_vector_t2>(); 
  expect_binary_match_return_t<F, vector<result_row_vector_t>, 
                               vector<row_vector_t1>, 
                               vector<row_vector_t2> >(); 
}

#endif
