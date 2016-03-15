#ifndef TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_REV_MAT_VECTORIZE_EXPECT_REV_STD_VECTOR_VALUE_HPP

#include <stan/math/rev/core/var.hpp>
#include <vector>
#include <test/unit/math/rev/mat/vectorize/build_rev_vector.hpp>
#include <test/unit/math/rev/mat/vectorize/expect_val_deriv_eq.hpp>

template <typename F>
void expect_rev_std_vector_value() {
  using std::vector;
  using stan::math::var;

  for (size_t i = 0; i < F::valid_inputs().size(); ++i) {
    vector<var> y = build_rev_vector<F>();
    vector<var> z = build_rev_vector<F>();
    vector<var> fz = F::template apply<vector<var> >(z);
    EXPECT_EQ(z.size(), fz.size());
    expect_val_deriv_eq(F::apply_base(y[i]), y[i], fz[i], z[i]);
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < F::valid_inputs().size(); ++j) {
      vector<vector<var> > a;
      vector<vector<var> > b;
      for (size_t i = 0; i < vector_vector_size; ++i) {
        a.push_back(build_rev_vector<F>());
        b.push_back(build_rev_vector<F>());
      } 
      vector<vector<var> > fb = F::template apply<vector<vector<var> > >(b);
      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      expect_val_deriv_eq(F::apply_base(a[i][j]), a[i][j], 
                          fb[i][j], b[i][j]);
    }
  }
}    

#endif
