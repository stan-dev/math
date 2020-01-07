#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, value_of_rec) {
  using stan::math::value_of_rec;
  using stan::math::var;

  var v_a(5.0);

  EXPECT_FLOAT_EQ(5.0, value_of_rec(v_a));
}
#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixRevArr, value_of_rec) {
  using stan::math::value_of_rec;
  using stan::math::var;
  using std::vector;

  vector<double> a_vals;
  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;
  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  vector<var> a;
  a = stan::math::to_var(a_vals);
  vector<var> b;
  b = stan::math::to_var(b_vals);

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val(), d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val(), d_a[i]);
}
