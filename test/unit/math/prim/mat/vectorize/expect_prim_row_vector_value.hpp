#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_ROW_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_ROW_VECTOR_VALUE_HPP

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <vector>
#include <stan/math/prim/scal/fun/is_nan.hpp>

template <typename F>
void expect_prim_row_vector_value() {
  using Eigen::RowVectorXd;
  using std::vector;
  using stan::math::is_nan;

  std::vector<double> valid_inputs = F::valid_inputs();

  RowVectorXd c = 
    RowVectorXd::Map(valid_inputs.data(), valid_inputs.size());

  RowVectorXd fc = F::template apply<RowVectorXd>(c);
  EXPECT_EQ(c.size(), fc.size());
  for (int i = 0; i < fc.size(); ++i) {
    if (is_nan(F::apply_base(c(i))) && is_nan(fc(i))) continue;
    EXPECT_FLOAT_EQ(F::apply_base(c(i)), fc(i));
  }

  vector<RowVectorXd> d;
  d.push_back(c);
  d.push_back(c);

  vector<RowVectorXd> fd = F::template apply<vector<RowVectorXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j) {
      if (is_nan(F::apply_base(d[i](j))) && is_nan(fd[i](j))) continue;
      EXPECT_FLOAT_EQ(F::apply_base(d[i](j)), fd[i](j));
    }
  }
}

#endif
