#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_eq.hpp>

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

    MatrixXvar b = build_matrix<F>(template_matrix);
    expect_eq(F::apply_base(b(i)), b(i), fa(i), a(i));
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_matrix.size(); ++j) {

      vector<MatrixXvar> c;
      for (size_t k = 0; k < vector_matrix_size; ++k)
        c.push_back(build_matrix<F>(template_matrix));
      vector<MatrixXvar> fc = F::template apply<vector<MatrixXvar> >(c);

      EXPECT_EQ(c[i].size(), fc[i].size());
      EXPECT_EQ(c[i].rows(), fc[i].rows());
      EXPECT_EQ(c[i].cols(), fc[i].cols());

      vector<MatrixXvar> d;
      for (size_t k = 0; k < vector_matrix_size; ++k)
        d.push_back(build_matrix<F>(template_matrix));
      expect_eq(F::apply_base(d[i](j)), d[i](j), fc[i](j), c[i](j));
    }
  }

  MatrixXvar e = build_matrix<F>(template_matrix);
  MatrixXvar feb = F::template apply<MatrixXvar>(e.block(1, 1, 1, 1));
  MatrixXvar f = build_matrix<F>(template_matrix).block(1, 1, 1, 1);
  expect_eq(F::apply_base(f(0,0)), f(0,0), feb(0,0), e(1,1));
}

#endif
