#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_ROW_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_ROW_VECTOR_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F>
void expect_prim_row_vector_value() {
  using Eigen::RowVectorXd;
  using std::vector;

  std::vector<double> valid_inputs = F::valid_inputs();

  RowVectorXd c = RowVectorXd::Map(valid_inputs.data(), valid_inputs.size());

  RowVectorXd fc = F::template apply<RowVectorXd>(c);
  EXPECT_EQ(c.size(), fc.size());
  for (int i = 0; i < fc.size(); ++i) {
    expect_val_eq(F::apply_base(c(i)), fc(i));
  }

  vector<RowVectorXd> d;
  d.push_back(c);
  d.push_back(c);

  vector<RowVectorXd> fd = F::template apply<vector<RowVectorXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j) {
      expect_val_eq(F::apply_base(d[i](j)), fd[i](j));
    }
  }
}

#endif
