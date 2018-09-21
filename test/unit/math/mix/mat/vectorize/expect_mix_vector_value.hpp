#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_VECTOR_VALUE_HPP

#include <test/unit/math/mix/mat/vectorize/build_mix_matrix.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_val_deriv_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename T>
void expect_mix_vector_value() {
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> vector_t;

  size_t num_inputs = F::valid_inputs().size();
  vector_t template_v(num_inputs);

  for (size_t i = 0; i < num_inputs; ++i) {
    vector_t a = build_mix_matrix<F>(template_v, i);
    vector_t b = build_mix_matrix<F>(template_v, i);
    vector_t fb = F::template apply<vector_t>(b);
    EXPECT_EQ(b.size(), fb.size());
    expect_val_deriv_eq(F::apply_base(a(i)), a(i), fb(i), b(i));
  }

  size_t vector_vector_size = 2;
  for (size_t i = 0; i < vector_vector_size; ++i) {
    for (size_t j = 0; j < num_inputs; ++j) {
      vector<vector_t> c;
      vector<vector_t> d;
      for (size_t k = 0; k < vector_vector_size; ++k)
        if (k == i) {
          c.push_back(build_mix_matrix<F>(template_v, j));
          d.push_back(build_mix_matrix<F>(template_v, j));
        } else {
          c.push_back(build_mix_matrix<F>(template_v));
          d.push_back(build_mix_matrix<F>(template_v));
        }
      vector<vector_t> fd = F::template apply<vector<vector_t> >(d);
      EXPECT_EQ(d.size(), fd.size());
      EXPECT_EQ(d[i].size(), fd[i].size());
      expect_val_deriv_eq(F::apply_base(c[i](j)), c[i](j), fd[i](j), d[i](j));
    }
  }
}

#endif
