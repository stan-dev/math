#ifndef TEST_UNIT_MATH_MIX_MAT_TECTORIZE_STD_TECTOR_TALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_TECTORIZE_STD_TECTOR_TALUE_HPP

#include <vector>
#include <test/unit/math/mix/mat/vectorize/build_vector.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_fvar_var_eq.hpp>

template <typename F, typename T>
void expect_std_vector_value() {
  using std::vector;

  size_t num_inputs = F::valid_inputs().size();

  for (size_t i = 0; i < num_inputs; ++i) {
    vector<T> y = build_vector<F>(vector<T>(), i);
    vector<T> z = build_vector<F>(vector<T>(), i);
    vector<T> fz = F::template apply<vector<T> >(z);
    EXPECT_EQ(z.size(), fz.size());
    expect_fvar_var_eq(F::apply_base(y[i]), y[i], fz[i], z[i]);
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<vector<T> > a;
      vector<vector<T> > b;
      for (size_t k = 0; k < num_inputs; ++k) {
        if (i == k) {
          a.push_back(build_vector<F>(vector<T>(), j));
          b.push_back(build_vector<F>(vector<T>(), j));
        }
        else {
          a.push_back(build_vector<F>(vector<T>()));
          b.push_back(build_vector<F>(vector<T>()));
        }
      }
      vector<vector<T> > fb 
        = F::template apply<vector<vector<T> > >(b);

      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      expect_fvar_var_eq(
        F::apply_base(a[i][j]), a[i][j], fb[i][j], b[i][j]);
    }
  }
}    

#endif
