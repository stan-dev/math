#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_ROW_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_ROW_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_rev_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_val_deriv_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F>
void expect_rev_row_vector_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, 1, Eigen::Dynamic> RowVectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  RowVectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    RowVectorXvar b = build_rev_matrix<F>(template_vector);
    RowVectorXvar c = build_rev_matrix<F>(template_vector);
    RowVectorXvar fc = F::template apply<RowVectorXvar>(c);
    EXPECT_EQ(c.size(), fc.size());
    expect_val_deriv_eq(F::apply_base(b(i)), b(i), fc(i), c(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<RowVectorXvar> d;
      vector<RowVectorXvar> e;
      for (size_t k = 0; k < num_inputs; ++k) {
        d.push_back(build_rev_matrix<F>(template_vector));
        e.push_back(build_rev_matrix<F>(template_vector));
      }
      vector<RowVectorXvar> fe = F::template apply<vector<RowVectorXvar> >(e);
      EXPECT_EQ(e[i].size(), fe[i].size());
      EXPECT_EQ(e.size(), fe.size());
      expect_val_deriv_eq(F::apply_base(d[i](j)), d[i](j), fe[i](j), e[i](j));
    }
  }
}

#endif
