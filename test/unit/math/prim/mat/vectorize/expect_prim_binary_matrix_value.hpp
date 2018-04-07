#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_BINARY_MATRIX_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <test/unit/math/prim/mat/vectorize/build_template_matrix.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F, typename matrix_t>
void expect_prim_binary_matrix_matrix_eq(const matrix_t& input_m1,
                                         const matrix_t& input_m2) {
  matrix_t fa = F::template apply<matrix_t>(input_m1, input_m2);
  EXPECT_EQ(input_m1.size(), fa.size());
  EXPECT_EQ(input_m2.size(), fa.size());
  for (int i = 0; i < input_m1.size(); ++i) {
    double exp_v = F::apply_base(input_m1(i), input_m2(i));
    expect_val_eq(exp_v, fa(i));
  }
}

template <typename F, typename matrix_t> void
expect_prim_binary_std_vector_matrix_std_vector_matrix_eq(
    const std::vector<matrix_t>& input_mv1,
    const std::vector<matrix_t>& input_mv2) {
  using std::vector;

  std::vector<matrix_t> fa = F::template apply<vector<matrix_t> >(
      input_mv1, input_mv2);
  EXPECT_EQ(input_mv1.size(), fa.size());
  EXPECT_EQ(input_mv2.size(), fa.size());
  for (size_t i = 0; i < input_mv1.size(); ++i) {
    EXPECT_EQ(input_mv1[i].size(), fa[i].size());
    EXPECT_EQ(input_mv2[i].size(), fa[i].size());
    for (int j = 0; j < input_mv1[i].size(); ++j) {
      double exp_v = F::apply_base(input_mv1[i](j), input_mv2[i](j));
      expect_val_eq(exp_v, fa[i](j));
    }
  }
}

template <typename F, typename matrix_t>
void expect_prim_binary_matrix_value(matrix_t template_m) {
  using Eigen::MatrixXd;
  using std::vector;

  vector<double> valid_inputs1 = F::valid_inputs1();
  vector<double> valid_inputs2 = F::valid_inputs2();

  matrix_t a1 = build_template_matrix(template_m, valid_inputs1.size(), 3);
  for (int i = 0; i < a1.size(); ++i) {
    a1(i) =  valid_inputs1[(i % valid_inputs1.size())];
  }

  matrix_t a2 = build_template_matrix(template_m, valid_inputs1.size(), 3);
  for (int i = 0; i < a2.size(); ++i) {
    a2(i) =  valid_inputs2[(i % valid_inputs2.size())];
  }

  expect_prim_binary_matrix_matrix_eq<F>(a1, a2);

  if (a1.rows() > 1 && a1.cols() > 1) {
    matrix_t ab1 = a1.block(1, 1, 1, 1);
    matrix_t ab2 = a2.block(1, 1, 1, 1);
    expect_prim_binary_matrix_matrix_eq<F>(ab1, ab2);
  }

  vector<MatrixXd> d1;
  d1.push_back(a1);
  d1.push_back(a1);
  vector<MatrixXd> d2;
  d2.push_back(a2);
  d2.push_back(a2);
  expect_prim_binary_std_vector_matrix_std_vector_matrix_eq<F>(d1, d2);
}

#endif
