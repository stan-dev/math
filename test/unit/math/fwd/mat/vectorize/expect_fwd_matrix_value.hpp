#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_FWD_MATRIX_VALUE_HPP

#include <stan/math/fwd/mat.hpp>
#include <test/unit/math/fwd/mat/vectorize/build_fwd_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_val_deriv_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename T>
void expect_fwd_matrix_value() {
  using stan::math::fvar;
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

  int num_inputs = F::valid_inputs().size();
  int num_cols = 3;
  matrix_t template_m(num_inputs, num_cols);

  for (int i = 0; i < template_m.size(); ++i) {
    matrix_t a = build_fwd_matrix<F>(template_m, i);
    matrix_t fa = F::template apply<matrix_t>(a);
    EXPECT_EQ(a.size(), fa.size());
    expect_val_deriv_eq(F::apply_base(a(i)), fa(i));
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_m.size(); ++j) {
      vector<matrix_t> b;
      for (size_t k = 0; k < vector_matrix_size; ++k)
        if (k == i)
          b.push_back(build_fwd_matrix<F>(template_m, j));
        else
          b.push_back(build_fwd_matrix<F>(template_m));
      vector<matrix_t> fb = F::template apply<vector<matrix_t> >(b);
      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      expect_val_deriv_eq(F::apply_base(b[i](j)), fb[i](j));
    }
  }

  int seed_i = num_inputs + 1;
  matrix_t a = build_fwd_matrix<F>(template_m, seed_i);
  matrix_t fab = F::template apply<matrix_t>(a.block(1, 1, 1, 1));
  expect_val_deriv_eq(F::apply_base(a(1, 1)), fab(0, 0));
}

#endif
