#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_VECTOR_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <Eigen/Dense>
#include <vector>

template <typename F>
void expect_prim_vector_value() {
  using Eigen::VectorXd;
  using std::vector;

  std::vector<double> valid_inputs = F::valid_inputs();

  VectorXd b = VectorXd::Map(valid_inputs.data(), valid_inputs.size());
  VectorXd fb = F::template apply<VectorXd>(b);
  EXPECT_EQ(b.size(), fb.size());
  for (int i = 0; i < fb.size(); ++i) {
    expect_val_eq(F::apply_base(b(i)), fb(i));
  }

  vector<VectorXd> d;
  d.push_back(b);
  d.push_back(b);
  vector<VectorXd> fd = F::template apply<vector<VectorXd> >(d);
  EXPECT_EQ(d.size(), fd.size());
  for (size_t i = 0; i < fd.size(); ++i) {
    EXPECT_EQ(d[i].size(), fd[i].size());
    for (int j = 0; j < fd[i].size(); ++j) {
      expect_val_eq(F::apply_base(d[i](j)), fd[i](j));
    }
  }
}
#endif
