#ifndef TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_MIX_MAT_VECTORIZE_EXPECT_MATRIX_VALUE_HPP

#include <vector>
#include <Eigen/Dense>
#include <test/unit/math/mix/mat/vectorize/build_matrix.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_eq.hpp>

template <typename F, typename T>
void expect_matrix_value() {
  using std::vector;
  typedef 
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

  size_t num_inputs = F::valid_inputs().size();
  size_t num_cols = 3;
  matrix_t template_matrix(num_inputs, num_cols);

  for (int i = 0; i < template_matrix.size(); ++i) {
    matrix_t y = build_matrix<F>(template_matrix, i);
    matrix_t z = build_matrix<F>(template_matrix, i);
    matrix_t fz = F::template apply<matrix_t>(z);
    EXPECT_EQ(z.size(), fz.size());
    expect_eq(F::apply_base(y(i)), y(i), fz(i), z(i));
  }

  size_t vector_matrix_size = 2;
  for (size_t i = 0; i < vector_matrix_size; ++i) {
    for (int j = 0; j < template_matrix.size(); ++j) {
      vector<matrix_t> a;
      vector<matrix_t> b;
      for (size_t k = 0; k < vector_matrix_size; ++k) {
        if (k == i) {
          a.push_back(build_matrix<F>(template_matrix, j));
          b.push_back(build_matrix<F>(template_matrix, j));
        } else {
          a.push_back(build_matrix<F>(template_matrix));
          b.push_back(build_matrix<F>(template_matrix));
        }
      }
      vector<matrix_t> fb = F::template apply<vector<matrix_t> >(b);
      EXPECT_EQ(b.size(), fb.size());
      EXPECT_EQ(b[i].size(), fb[i].size());
      expect_eq(
        F::apply_base(a[i](j)), a[i](j), fb[i](j), b[i](j));
    }
  }

  int block_i = 1;
  int block_j = 1;
  int seed_i = block_j * num_inputs + block_i; 
  matrix_t c = build_matrix<F>(template_matrix, seed_i);
  matrix_t d = build_matrix<F>(template_matrix, seed_i);
  matrix_t fab = foo(d.block(block_i, block_j, 1, 1));
  expect_eq(
    F::apply_base(c(block_i, block_j)), c(block_i, block_j),
              fab(0,0), d(block_i, block_j));
}

#endif
