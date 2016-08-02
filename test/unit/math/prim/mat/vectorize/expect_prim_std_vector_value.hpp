#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_STD_VECTOR_VALUE_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_PRIM_STD_VECTOR_VALUE_HPP

#include <test/unit/math/prim/mat/vectorize/expect_val_eq.hpp>
#include <vector>

template <typename F>
void expect_prim_std_vector_value() {
  using std::vector;

  vector<double> valid_inputs = F::valid_inputs();

  vector<double> fy = F::template apply<vector<double> >(valid_inputs);
  EXPECT_EQ(valid_inputs.size(), fy.size());
  for (size_t i = 0; i < valid_inputs.size(); ++i) {
    expect_val_eq(F::apply_base(valid_inputs[i]), fy[i]);
  }

  vector<vector<double> > z;
  z.push_back(valid_inputs);
  z.push_back(valid_inputs);
  vector<vector<double> > fz
    = F::template apply<vector<vector<double> > >(z);
  EXPECT_EQ(z.size(), fz.size());
  for (size_t i = 0; i < fz.size(); ++i) {
    EXPECT_EQ(z[i].size(), fz[i].size());
    for (size_t j = 0; j < z[i].size(); ++j) {
      expect_val_eq(F::apply_base(z[i][j]), fz[i][j]);
    }
  }

  vector<int> int_valid_inputs = F::int_valid_inputs();
  vector<double> u2 = F::template apply<vector<double> >(int_valid_inputs);
  EXPECT_EQ(int_valid_inputs.size(), u2.size());
  for (size_t i = 0; i < int_valid_inputs.size(); ++i) {
    expect_val_eq(F::apply_base(int_valid_inputs[i]), u2[i]);
  }
}

#endif
