#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATRIX_ERROR_HPP

#include <test/unit/math/prim/mat/vectorize/build_template_matrix.hpp>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <stdexcept>
#include <vector>

template <typename F, typename matrix_t, typename matrix_d>
void expect_binary_matrix_size_error(matrix_t template_tm,
                                     matrix_d template_dm) {
  using std::vector;

  matrix_t badsize_tm1 = build_template_matrix(template_tm, 4, 5);
  matrix_t badsize_tm2 = build_template_matrix(template_tm, 7, 5);
  matrix_t badsize_tm3 = build_template_matrix(template_tm, 4, 9);
  matrix_t badsize_tm4 = build_template_matrix(template_tm, 3, 6);
  matrix_t badsize_tm5 = build_template_matrix(template_tm, 5, 4);
  matrix_d badsize_dm1 = build_template_matrix(template_dm, 6, 13);
  matrix_d badsize_dm2 = build_template_matrix(template_dm, 4, 8);
  matrix_d badsize_dm3 = build_template_matrix(template_dm, 12, 5);
  matrix_d badsize_dm4 = build_template_matrix(template_dm, 5, 4);

  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm1),
               std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_dm1, badsize_tm1),
               std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm2),
               std::invalid_argument);
  EXPECT_THROW(F::template apply<matrix_t>(badsize_tm2, badsize_tm1),
               std::invalid_argument);

  if (badsize_tm1.rows() > 1 && badsize_tm1.cols() > 1) {
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm2),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_dm2, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm3),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_dm3, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_dm4),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_dm4, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm3),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm3, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm4),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm4, badsize_tm1),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm1, badsize_tm5),
                 std::invalid_argument);
    EXPECT_THROW(F::template apply<matrix_t>(badsize_tm5, badsize_tm1),
                 std::invalid_argument);
  }
}

template <typename F, typename result_t, typename T1, 
          typename T2, int R, int C>
void expect_binary_matrix_err_throw(
  Eigen::Matrix<T1, R, C> template_m1, 
  Eigen::Matrix<T2, R, C> template_m2,
  std::vector<double> invalid_inputs1,
  std::vector<double> invalid_inputs2) {

  typedef typename Eigen::Matrix<T1, R, C> matrix_t1;
  typedef typename Eigen::Matrix<T2, R, C> matrix_t2;
  typedef typename Eigen::Matrix<result_t, R, C> result_tm;

  for (size_t i = 0; i < invalid_inputs1.size(); ++i) {
    matrix_t1 a = Eigen::Matrix<T1, R, C>::Constant(
                    template_m1.rows(), template_m1.cols(), 
                    invalid_inputs1[i]); 
    matrix_t2 b = Eigen::Matrix<T2, R, C>::Constant(
                    template_m2.rows(), template_m2.cols(), 
                    invalid_inputs2[i]); 
    EXPECT_THROW(F::template apply<result_tm>(a, b), std::domain_error);

    matrix_t2 c = Eigen::Matrix<T2, R, C>::Constant(
                    template_m2.rows(), template_m2.cols(), 
                    invalid_inputs1[i]); 
    matrix_t1 d = Eigen::Matrix<T1, R, C>::Constant(
                    template_m1.rows(), template_m1.cols(), 
                    invalid_inputs2[i]); 
    EXPECT_THROW(F::template apply<result_tm>(c, d), std::domain_error);

    if (template_m1.rows() > 1 && template_m1.cols() > 1) {
      EXPECT_THROW(F::template apply<result_tm>(a.block(1, 1, 1, 1),
                                                b.block(1, 1, 1, 1)),
                   std::domain_error);
      EXPECT_THROW(F::template apply<result_tm>(c.block(1, 1, 1, 1),
                                                d.block(1, 1, 1, 1)),
                   std::domain_error);
    }
  }
}

template <typename F, typename result_t, typename T1, 
          typename T2, int R, int C>
void expect_binary_std_vector_matrix_err_throw(
  Eigen::Matrix<T1, R, C> template_m1, 
  Eigen::Matrix<T2, R, C> template_m2,
  std::vector<double> invalid_inputs1,
  std::vector<double> invalid_inputs2) {

  using std::vector;
  typedef typename Eigen::Matrix<T1, R, C> matrix_t1;
  typedef typename Eigen::Matrix<T2, R, C> matrix_t2;
  typedef typename Eigen::Matrix<result_t, R, C> result_tm;

  vector<matrix_t1> e;
  vector<matrix_t2> f;
  vector<matrix_t2> g;
  vector<matrix_t1> h;
  for (size_t i = 0; i < invalid_inputs1.size(); ++i) {
    e.clear();
    f.clear();
    g.clear();
    h.clear();
    matrix_t1 a = Eigen::Matrix<T1, R, C>::Constant(
                    template_m1.rows(), template_m1.cols(), 
                    invalid_inputs1[i]); 
    matrix_t2 b = Eigen::Matrix<T2, R, C>::Constant(
                    template_m2.rows(), template_m2.cols(), 
                    invalid_inputs2[i]); 
    matrix_t2 c = Eigen::Matrix<T2, R, C>::Constant(
                    template_m2.rows(), template_m2.cols(), 
                    invalid_inputs1[i]); 
    matrix_t1 d = Eigen::Matrix<T1, R, C>::Constant(
                    template_m1.rows(), template_m1.cols(), 
                    invalid_inputs2[i]); 
    for (size_t j = 0; j < 2; ++j) {
      e.push_back(a);
      f.push_back(b);
      g.push_back(c);
      h.push_back(d);
    }
    EXPECT_THROW(F::template apply<vector<result_tm> >(e, f), 
                 std::domain_error);
    EXPECT_THROW(F::template apply<vector<result_tm> >(g, h), 
                 std::domain_error);
  }
}

template <typename F, typename T, int R, int C>
void expect_binary_matrix_value_error(
  Eigen::Matrix<T, R, C> template_tm, 
  Eigen::Matrix<double, R, C> template_dm) {

  using std::vector;

  vector<double> invalid_inputs1 = F::invalid_inputs1();
  if (invalid_inputs1.size() == 0) return;
  vector<double> invalid_inputs2 = F::invalid_inputs2();

  template_dm = build_template_matrix(template_dm, 5, 3);
  template_tm = build_template_matrix(template_tm, 5, 3);

  expect_binary_matrix_err_throw<F, T>(template_tm, template_dm,
                                              invalid_inputs1,
                                              invalid_inputs2);
  expect_binary_matrix_err_throw<F, T>(template_tm, template_tm,
                                              invalid_inputs1,
                                              invalid_inputs2);

  expect_binary_std_vector_matrix_err_throw<F, T>(
    template_tm, template_dm, invalid_inputs1, invalid_inputs2);
  expect_binary_std_vector_matrix_err_throw<F, T>(
    template_tm, template_tm, invalid_inputs1, invalid_inputs2);
}


template <typename F, typename T, int R, int C>
void expect_binary_matrix_error(Eigen::Matrix<double, R, C> template_m) {
  Eigen::Matrix<T, R, C> template_tm;
  expect_binary_matrix_size_error<F>(template_tm, template_m);
  expect_binary_matrix_value_error<F>(template_tm, template_m);
}

#endif
