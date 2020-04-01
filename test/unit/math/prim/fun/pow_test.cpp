#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

namespace test1 {
  using namespace stan::math;
  void test_pow() {
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }
}

namespace test2 {
  using namespace std;
  using namespace stan::math;
  void test_pow() {
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }
}

namespace test3 {
  using namespace stan::math;
  void test_pow() {
    using std::pow;
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }
}

namespace test4 {
  using namespace std;
  void test_pow() {
    using stan::math::pow;
    EXPECT_FLOAT_EQ(pow(1.5, 1.5), 1.837117);
    EXPECT_FLOAT_EQ(pow(1.7, 2), 2.89);
  }
}

namespace test5 {
  using namespace std;
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
    using std::pow;
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
  test3::test_pow();
  test4::test_pow();
  test5::test_pow();
}
