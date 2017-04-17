#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_rev_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_val_deriv_eq.hpp>

template <typename F>
void expect_rev_vector_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, Eigen::Dynamic, 1> VectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  VectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    VectorXvar b = build_rev_matrix<F>(template_vector);
    VectorXvar c = build_rev_matrix<F>(template_vector);
    VectorXvar fc = F::template apply<VectorXvar>(c);
    EXPECT_EQ(c.size(), fc.size());
    expect_val_deriv_eq(F::apply_base(b(i)), b(i), fc(i), c(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<VectorXvar> d;
      vector<VectorXvar> e;
      for (size_t k = 0; k < num_inputs; ++k) {
        d.push_back(build_rev_matrix<F>(template_vector));
        e.push_back(build_rev_matrix<F>(template_vector));
      }
      vector<VectorXvar> fe = F::template apply<vector<VectorXvar> >(e);
      EXPECT_EQ(e[i].size(), fe[i].size());
      EXPECT_EQ(e.size(), fe.size());
      expect_val_deriv_eq(F::apply_base(d[i](j)), d[i](j), 
                          fe[i](j), e[i](j)); 
    }
  }
}

#endif
