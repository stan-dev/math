#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_MATRIX_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_rev_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_val_deriv_eq.hpp>

template <typename F>
void expect_rev_matrix_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> MatrixXvar;
  
  size_t num_cols = 3;
  size_t num_inputs = F::valid_inputs().size();
  MatrixXvar template_m(num_inputs, num_cols);

  for (int i = 0; i < template_m.size(); ++i) {
    MatrixXvar a = build_rev_matrix<F>(template_m);
    MatrixXvar b = build_rev_matrix<F>(template_m);
    MatrixXvar fb = F::template apply<MatrixXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    expect_val_deriv_eq(F::apply_base(a(i)), a(i), fb(i), b(i));
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_m.size(); ++j) {
      vector<MatrixXvar> c;
      vector<MatrixXvar> d;
      for (size_t k = 0; k < vector_matrix_size; ++k) {
        c.push_back(build_rev_matrix<F>(template_m));
        d.push_back(build_rev_matrix<F>(template_m));
      }
      vector<MatrixXvar> fd = F::template apply<vector<MatrixXvar> >(d);
      EXPECT_EQ(d[i].size(), fd[i].size());
      EXPECT_EQ(d[i].rows(), fd[i].rows());
      EXPECT_EQ(d[i].cols(), fd[i].cols());
      expect_val_deriv_eq(F::apply_base(c[i](j)), c[i](j),
                          fd[i](j), d[i](j));
    }
  }

  MatrixXvar e = build_rev_matrix<F>(template_m).block(1, 1, 1, 1);
  MatrixXvar f = build_rev_matrix<F>(template_m);
  MatrixXvar ffb = F::template apply<MatrixXvar>(f.block(1, 1, 1, 1));
  expect_val_deriv_eq(F::apply_base(e(0,0)), e(0,0), ffb(0,0), f(1,1));
}

#endif
