#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

// Make sure stan::math has a pow implementation
namespace test1 {
  using namespace stan::math;
  void test_pow() {
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }
}

// Make sure std and stan::math don't conflict
namespace test2 {
  using namespace std;
  using namespace stan::math;
  void test_pow() {
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }
}

TEST(MathFunctions, powNamespaceChecks) {
  EXPECT_FLOAT_EQ(stan::math::pow(1.5, 1.5), 1.837117);
  EXPECT_FLOAT_EQ(stan::math::pow(1.7, 2), 2.89);

  {
    using stan::math::pow;
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }

  {
    using stan::math::pow;
    using std::pow;
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }

  test1::test_pow();
  test2::test_pow();
}
