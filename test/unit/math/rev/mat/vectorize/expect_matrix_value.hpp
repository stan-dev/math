#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_autodiff.hpp>

template <typename F>
void expect_matrix_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> MatrixXvar;
  
  size_t num_cols = 3;
  size_t num_inputs = F::valid_inputs().size();
  MatrixXvar template_matrix(num_inputs, num_cols);

  for (int i = 0; i < template_matrix.size(); ++i) {
    MatrixXvar a = build_matrix<F>(template_matrix);
    MatrixXvar fa = F::template apply<MatrixXvar>(a);
    EXPECT_EQ(a.size(), fa.size());
    EXPECT_FLOAT_EQ(F::apply_base(a(i)).val(), fa(i).val());

    fa(i).grad();
    expect_autodiff<F>(a(i).val(), a(i).adj());
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_matrix.size(); ++j) {

      vector<MatrixXvar> b;
      for (size_t k = 0; k < vector_matrix_size; ++k)
        b.push_back(build_matrix<F>(template_matrix));
      vector<MatrixXvar> fb = F::template apply<vector<MatrixXvar> >(b);

      EXPECT_EQ(b[i].size(), fb[i].size());
      EXPECT_EQ(b[i].rows(), fb[i].rows());
      EXPECT_EQ(b[i].cols(), fb[i].cols());
      EXPECT_FLOAT_EQ(F::apply_base(b[i](j)).val(), fb[i](j).val()); 

      fb[i](j).grad();
      expect_autodiff<F>(b[i](j).val(), b[i](j).adj());
    }
  }

  MatrixXvar a = build_matrix<F>(template_matrix);
  MatrixXvar fab = F::template apply<MatrixXvar>(a.block(1, 1, 1, 1));
  EXPECT_FLOAT_EQ(F::apply_base(a(1,1)).val(), fab(0,0).val());
  fab(0,0).grad();
  expect_autodiff<F>(a(1,1).val(), a(1,1).adj());
}

#endif
