#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_VECTOR_VALUE_HPP

#include <stan/math/fwd/mat.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename T>
void expect_fwd_vector_value() {
  using stan::math::fvar;
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> vector_t;

  size_t num_inputs = F::valid_inputs().size();
  vector_t template_v(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    vector_t b = build_fwd_matrix<F>(template_v, i);
    vector_t fb = F::template apply<vector_t>(b);
    EXPECT_EQ(b.size(), fb.size());
    expect_val_deriv_eq(F::apply_base(b(i)), fb(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<vector_t> c;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i)
          c.push_back(build_fwd_matrix<F>(template_v, j));
        else
          c.push_back(build_fwd_matrix<F>(template_v));
      vector<vector_t> fc = F::template apply<vector<vector_t> >(c);
      EXPECT_EQ(c.size(), fc.size());
      EXPECT_EQ(c[i].size(), fc[i].size());
      expect_val_deriv_eq(F::apply_base(c[i](j)), fc[i](j));
    }
  }
}

#endif
