#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

// Several special cases can't be estimated by Boost's pFq for reference
// test values (as they fail convergence checks for the infinite sums) and are
// taken from WolframAlpha instead
TEST(MathFunctions, hypergeometric_2F1_special_cases) {
  using Eigen::VectorXd;
  using stan::math::hypergeometric_2F1;
  using stan::math::hypergeometric_pFq;
  using stan::math::inv;

  VectorXd a(2);
  VectorXd b(1);
  a << 1, 1;
  b << 2;
  double z = -5;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B1%2C1%2C2%2C-5%5D
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z), std::log(6) / 5);

  a << 4, 4;
  b << 10;
  z = 1;

  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z),
                  hypergeometric_pFq(a, b, z));

  a << 6, 4;
  b << 4;
  z = 8;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B6%2C4%2C4%2C8%5D
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z), inv(117649));

  a << 4, 6;
  b << 4;
  z = 8;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B6%2C4%2C4%2C8%5D
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z), inv(117649));

  a << 0.5, 0.5;
  b << 1.5;
  z = 0.8;
  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B1%2F2%2C1%2F2%2C3%2F2%2C0.8%5D
  // // NOLINT
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z),
                  1.2378298970946586875911981729041590167);

  a << 0.5, 0.5;
  b << 1.5;
  z = -0.8;
  // https://www.wolframalpha.com/input?i=Hypergeometric2F1[1%2F2%2C1%2F2%2C3%2F2%2C-0.8]
  // // NOLINT
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z),
                  0.8997031444420006140233288415980043158);

  a << 4, 4;
  b << 10;
  z = 0.0;

  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z),
                  hypergeometric_pFq(a, b, z));

  a << 1, 4;
  b << 4;
  z = -0.8;

  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z),
                  hypergeometric_pFq(a, b, z));

  a << 1.5, 2.0;
  b << 3.0;
  z = -1.5;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B3%2F2%2C+2%2C+3%2C+-1.5%5D
  // // NOLINT
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z),
                  0.3797233104317609464875119664051608419);

  a << 1.5, 2.0;
  b << 5.5;
  z = 1.0;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B1.5%2C+2.0%2C+5.5%2C+1.0%5D
  // // NOLINT
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z), 2.625);

  a << 4, 4;
  b << 5;
  z = 0.5;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B4%2C4%2C5%2C0.5%5D
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z),
                  8.9719137774968335306304775600100329765);

  a << 4, 6;
  b << 5;
  z = 4;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B4%2C6%2C5%2C4%5D
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z), -inv(1215));

  a << 4, 1;
  b << 0.5;
  z = -5;

  // https://www.wolframalpha.com/input?i=Hypergeometric2F1%5B4%2C6%2C5%2C4%5D
  EXPECT_FLOAT_EQ(hypergeometric_2F1(a[0], a[1], b[0], z), -0.0399473);
}
