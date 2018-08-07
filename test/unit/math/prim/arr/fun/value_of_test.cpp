#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <type_traits>

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

TEST(MathFunctions, value_of_return_type_short_circuit) {
  std::vector<double> a(5, 0);
  bool r1 = std::is_same<decltype(stan::math::value_of(a)),
                         std::vector<double>>::value;
  EXPECT_FALSE(r1);
  bool r2 = std::is_same<decltype(stan::math::value_of(a)),
                         std::vector<double>&>::value;
  EXPECT_FALSE(r2);
  bool r3 = std::is_same<decltype(stan::math::value_of(a)),
                         const std::vector<double>>::value;
  EXPECT_FALSE(r3);
  bool r4 = std::is_same<decltype(stan::math::value_of(a)),
                         const std::vector<double>&>::value;
  EXPECT_TRUE(r4);
}
