#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_STD_VECTOR_VALUE_HPP

#include <stan/math/fwd/mat.hpp>
#include <vector>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_vector.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>

template <typename F, typename T>
void expect_fwd_std_vector_value() {
  using std::vector;

  size_t num_inputs = F::valid_inputs().size();
  for (size_t i = 0; i < num_inputs; ++i) {
    vector<T> y = build_fwd_vector<F>(vector<T>(), i);
    vector<T> fy = F::template apply<vector<T> >(y);
    EXPECT_EQ(y.size(), fy.size());
    expect_val_deriv_eq(F::apply_base(y[i]), fy[i]);
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<vector<T> > z;
      for (size_t k = 0; k < num_inputs; ++k) {
        if (i == k)
          z.push_back(build_fwd_vector<F>(vector<T>(), j));
        else
          z.push_back(build_fwd_vector<F>(vector<T>()));
      }
      vector<vector<T> > fz = 
        F::template apply<vector<vector<T> > >(z);
      EXPECT_EQ(z.size(), fz.size());
      EXPECT_EQ(z[i].size(), fz[i].size());
      expect_val_deriv_eq(F::apply_base(z[i][j]), fz[i][j]);
    }
  }
}    

#endif
