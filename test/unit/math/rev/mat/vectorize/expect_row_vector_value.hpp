#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_ROW_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_ROW_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_eq.hpp>

template <typename F>
void expect_row_vector_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, 1, Eigen::Dynamic> RowVectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  RowVectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    RowVectorXvar b = build_matrix<F>(template_vector);
    RowVectorXvar fb = F::template apply<RowVectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    
    RowVectorXvar c = build_matrix<F>(template_vector);
    expect_eq(F::apply_base(c(i)), c(i), fb(i), b(i));
  }

  size_t vector_vector_size = 2;

  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<RowVectorXvar> d;
      for (size_t k = 0; k < num_inputs; ++k)
        d.push_back(build_matrix<F>(template_vector));
      vector<RowVectorXvar> fd = 
        F::template apply<vector<RowVectorXvar> >(d);
      EXPECT_EQ(d[i].size(), fd[i].size());
      EXPECT_EQ(d.size(), fd.size());

      vector<RowVectorXvar> e;
      for (size_t k = 0; k < num_inputs; ++k)
        e.push_back(build_matrix<F>(template_vector));
      expect_eq(F::apply_base(e[i](j)), e[i](j), fd[i](j), d[i](j)); 
    }
  }
}

#endif
