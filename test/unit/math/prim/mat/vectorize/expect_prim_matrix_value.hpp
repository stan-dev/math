#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_MATRIX_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_MATRIX_VALUE_HPP

#include <vector>
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/is_nan.hpp>

template <typename F>
void expect_prim_matrix_value() {
  using Eigen::MatrixXd;
  using std::vector;
  using stan::math::is_nan;

  vector<double> valid_inputs = F::valid_inputs();
  MatrixXd a(valid_inputs.size(), 3);

  for (int i = 0; i < a.size(); i++) {
    a(i) =  valid_inputs[(i % valid_inputs.size())];
  }

  MatrixXd fa = F::template apply<MatrixXd>(a);
  EXPECT_EQ(a.size(), fa.size());
  for (int i = 0; i < fa.size(); ++i) {
    if (is_nan(F::apply_base(a(i))) && is_nan(fa(i))) continue;
    EXPECT_FLOAT_EQ(F::apply_base(a(i)), fa(i));
  }

  vector<MatrixXd> d;
  d.push_back(a);
  d.push_back(a);
  vector<MatrixXd> fd = F::template apply<vector<MatrixXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].rows(), fd[i].rows());
    EXPECT_EQ(d[i].cols(), fd[i].cols());
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j) {
      if (is_nan(F::apply_base(d[i](j))) && is_nan(fd[i](j))) return;
      EXPECT_FLOAT_EQ(F::apply_base(d[i](j)), fd[i](j));
    }
  }

  MatrixXd fab = F::template apply<MatrixXd>(a.block(1, 1, 1, 1));
  if (is_nan(F::apply_base(a(1,1))) && is_nan(fab(0,0))) return;
  EXPECT_FLOAT_EQ(F::apply_base(a(1,1)), fab(0,0));
}

#endif
