#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <Eigen/Dense>
#include <vector>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_autodiff.hpp>

template <typename F>
void expect_vector_value() {
  using stan::math::var;
  using std::vector;
  typedef Eigen::Matrix<var, Eigen::Dynamic, 1> VectorXvar;

  size_t num_inputs = F::valid_inputs().size();
  VectorXvar template_vector(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    VectorXvar b = build_matrix<F>(template_vector);
    VectorXvar fb = F::template apply<VectorXvar>(b);
    EXPECT_EQ(b.size(), fb.size());
    
    EXPECT_FLOAT_EQ(F::apply_base(b(i)).val(), fb(i).val());

    fb(i).grad();
    expect_autodiff<F>(b(i).val(), b(i).adj());
  }

  size_t vector_vector_size = 2;

  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<VectorXvar> vb;
      for (size_t k = 0; k < num_inputs; ++k)
        vb.push_back(build_matrix<F>(template_vector));
      vector<VectorXvar> fvb = F::template apply<vector<VectorXvar> >(vb);
      EXPECT_EQ(vb[i].size(), fvb[i].size());
      EXPECT_EQ(vb.size(), fvb.size());
      EXPECT_FLOAT_EQ(F::apply_base(vb[i](j)).val(), fvb[i](j).val());

      fvb[i](j).grad();
      expect_autodiff<F>(vb[i](j).val(), vb[i](j).adj()); 
    }
  }
}

#endif
