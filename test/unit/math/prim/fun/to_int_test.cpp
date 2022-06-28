#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, to_int) {
  using stan::math::to_int;

  EXPECT_EQ(2, to_int(2.0));
  EXPECT_EQ(2, to_int(2.1));
  EXPECT_EQ(2, to_int(2.9));
  EXPECT_EQ(2, to_int(2.999999999));

  EXPECT_EQ(-36574, to_int(-36574.0));
  EXPECT_EQ(-36574, to_int(-36574.1));
  EXPECT_EQ(-36574, to_int(-36574.9));
  EXPECT_EQ(-36574, to_int(-36574.999999999));

  EXPECT_THROW(to_int(std::numeric_limits<int>::max() + 1.0),
               std::invalid_argument);
  EXPECT_THROW(to_int(std::numeric_limits<int>::min() - 1.0),
               std::invalid_argument);
}

TEST(MathFunctions, to_int_vec) {
  using stan::math::to_int;

  std::vector<double> inputs_1{2.1, -34.64, 10.89, 1000000};
  std::vector<double> inputs_2{-409831.987, 403.1, 10.61, -0.00001};
  std::vector<std::vector<double>> inputs{inputs_1, inputs_2};

  std::vector<int> target_result_1{2, -34, 10, 1000000};
  std::vector<int> target_result_2{-409831, 403, 10, 0};
  std::vector<std::vector<int>> target_result{target_result_1, target_result_2};

  EXPECT_STD_VECTOR_EQ(to_int(inputs), target_result);

  inputs[0][2] = std::numeric_limits<int>::min() - 1.0;
  EXPECT_THROW(to_int(inputs), std::invalid_argument);
}
