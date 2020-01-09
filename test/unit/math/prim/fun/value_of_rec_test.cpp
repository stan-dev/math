#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, value_of_rec) {
  using stan::math::value_of_rec;
  double x = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of_rec(x));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(5));
}
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixPrimArr, value_of_rec) {
  using stan::math::value_of_rec;
  using std::vector;

  vector<double> a;
  for (size_t i = 0; i < 10; ++i)
    a.push_back(i + 1);

  vector<double> b;
  for (size_t i = 10; i < 15; ++i)
    b.push_back(i + 1);

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i], d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i], d_a[i]);
}
