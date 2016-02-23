#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP

#include <stan/math/fwd/core/fvar.hpp>
#include <vector>
#include <Eigen/Dense>
#include <test/unit/math/fwd/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/fwd/mat/vectorize/expect_eq.hpp>

template <typename F, typename T>
void expect_matrix_value() {
  using stan::math::fvar;
  using std::vector;
  typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

  int num_inputs = F::valid_inputs().size();
  int num_cols = 3;
  matrix_t template_matrix(num_inputs, num_cols);

  for (int i = 0; i < template_matrix.size(); ++i) {
    matrix_t a = build_matrix<F>(template_matrix, i);
    matrix_t fa = F::template apply<matrix_t>(a);
    EXPECT_EQ(a.size(), fa.size());
    expect_eq(F::apply_base(a(i)), fa(i));
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_matrix.size(); ++j) {
      vector<matrix_t> b;
      for (size_t k = 0; k < vector_matrix_size; ++k)
        if (k == i)
          b.push_back(build_matrix<F>(template_matrix, j));
        else
          b.push_back(build_matrix<F>(template_matrix));
      vector<matrix_t> fb = F::template apply<vector<matrix_t> >(b);
      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      expect_eq(F::apply_base(b[i](j)), fb[i](j));
    }
  }

  int block_i = 1;
  int block_j = 1;
  int seed_i = block_j * num_inputs + block_i;
  matrix_t a = build_matrix<F>(template_matrix, seed_i);
  matrix_t fab = foo(a.block(block_i, block_j, 1, 1));
  expect_eq(F::apply_base(a(1,1)), fab(0,0));
}

#endif
