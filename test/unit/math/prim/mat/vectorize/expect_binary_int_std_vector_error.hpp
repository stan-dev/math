#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_INT_STD_VECTOR_ERROR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_INT_STD_VECTOR_ERROR_HPP

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

template <typename F>
void expect_binary_int_std_vector_error() {
  using std::vector;

  vector<int> a1(5);
  vector<int> a2(9);
  EXPECT_THROW(F::template apply<vector<double> >(a1, a2),
               std::invalid_argument);

  vector<int> invalid_inputs1 = F::int_invalid_inputs1();
  vector<int> invalid_inputs2 = F::int_invalid_inputs2();
  if (invalid_inputs1.size() == 0)
    return;
  EXPECT_THROW(
      F::template apply<vector<double> >(invalid_inputs1, invalid_inputs2),
      std::domain_error);

  vector<vector<int> > z1;
  z1.push_back(invalid_inputs1);
  z1.push_back(invalid_inputs1);
  vector<vector<int> > z2;
  z2.push_back(invalid_inputs2);
  z2.push_back(invalid_inputs2);
  EXPECT_THROW(F::template apply<vector<vector<double> > >(z1, z2),
               std::domain_error);
}

#endif
