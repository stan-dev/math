#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/arr/util.hpp>
#include <vector>

TEST(AgradRev, sum_std_vector) {
  using stan::math::sum;
  using stan::math::var;
  using std::vector;

  vector<var> x;
  for (size_t i = 0; i < 6; ++i)
    x.push_back(i + 1);

  var fx = 3.7 * sum(x);
  EXPECT_FLOAT_EQ(3.7 * 21.0, fx.val());

  vector<double> gx;
  fx.grad(x, gx);
  EXPECT_EQ(6, gx.size());
  for (size_t i = 0; i < 6; ++i)
    EXPECT_FLOAT_EQ(3.7, gx[i]);

  x = vector<var>();
  EXPECT_FLOAT_EQ(0.0, sum(x).val());
}

TEST(AgradRev, check_varis_on_stack) {
  std::vector<stan::math::var> x;
  for (size_t i = 0; i < 6; ++i)
    x.push_back(i + 1);
  test::check_varis_on_stack(stan::math::sum(x));
}
