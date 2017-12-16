#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MIX_MATRIX_VALUE_HPP

#include <test/unit/math/mix/mat/vectorize/build_mix_matrix.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_val_deriv_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename T>
void expect_mix_matrix_value() {
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

  size_t num_inputs = F::valid_inputs().size();
  size_t num_cols = 3;
  matrix_t template_m(num_inputs, num_cols);

  for (int i = 0; i < template_m.size(); ++i) {
    matrix_t y = build_mix_matrix<F>(template_m, i);
    matrix_t z = build_mix_matrix<F>(template_m, i);
    matrix_t fz = F::template apply<matrix_t>(z);
    EXPECT_EQ(z.size(), fz.size());
    expect_val_deriv_eq(F::apply_base(y(i)), y(i), fz(i), z(i));
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_m.size(); ++j) {
      vector<matrix_t> a;
      vector<matrix_t> b;
      for (size_t k = 0; k < vector_matrix_size; ++k) {
        if (k == i) {
          a.push_back(build_mix_matrix<F>(template_m, j));
          b.push_back(build_mix_matrix<F>(template_m, j));
        } else {
          a.push_back(build_mix_matrix<F>(template_m));
          b.push_back(build_mix_matrix<F>(template_m));
        }
      }
      vector<matrix_t> fb = F::template apply<vector<matrix_t> >(b);
      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      expect_val_deriv_eq(F::apply_base(a[i](j)), a[i](j), fb[i](j), b[i](j));
    }
  }

  int seed_i = num_inputs + 1;
  matrix_t c = build_mix_matrix<F>(template_m, seed_i);
  matrix_t d = build_mix_matrix<F>(template_m, seed_i);
  matrix_t fab = F::template apply<matrix_t>(d.block(1, 1, 1, 1));
  expect_val_deriv_eq(F::apply_base(c(1, 1)), c(1, 1), fab(0, 0), d(1, 1));
}

#endif
