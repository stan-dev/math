#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrix, value_of) {
  using stan::math::value_of;
  using std::vector;

  vector<double> a;
  for (size_t i = 0; i < 10; ++i)
    a.push_back(i + 1);

  vector<double> b;
  for (size_t i = 10; i < 15; ++i)
    b.push_back(i + 1);

  vector<double> d_a = value_of(a);
  vector<double> d_b = value_of(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i], d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i], d_a[i]);
}
