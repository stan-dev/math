#include <stan/math/prim/mat.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(MathFunctions, adapt_double) {
  double x = 1.0;

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(1, a.size());
  EXPECT_FLOAT_EQ(x, a(0));

  a(0) = 5.0;
  EXPECT_FLOAT_EQ(5.0, a(0));
}

TEST(MathFunctions, adapt_int) {
  int x = 1;

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(0, a.size());
}

TEST(MathFunctions, adapt_std_vector_double) {
  std::vector<double> x = {{1.0, 2.0}};

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(x.size(), a.size());
  for (size_t i = 0; i < x.size(); ++i)
    EXPECT_FLOAT_EQ(x[i], a(i));
}

TEST(MathFunctions, adapt_std_vector_int) {
  std::vector<int> x = {{1, 2}};

  auto a = stan::math::make_variable_adapter<double>(x);

  EXPECT_EQ(0, a.size());
}

TEST(MathFunctions, adapt_all_types) {
  double xd = 1.0;
  int xi = 1;
  std::vector<double> xdv = {{1.0, 2.0}};
  std::vector<int> xiv = {{1, 2}};

  auto a = stan::math::make_variable_adapter<double>(xd, xi, xdv, xiv, xd);

  EXPECT_EQ(2 + xdv.size(), a.size());
  EXPECT_FLOAT_EQ(xd, a(0));
  for (size_t i = 0; i < xdv.size(); ++i)
    EXPECT_FLOAT_EQ(xdv[i], a(1 + i));
  EXPECT_FLOAT_EQ(xd, a(1 + xdv.size()));

  a(xdv.size()) = 5.0;
  EXPECT_FLOAT_EQ(5.0, a(xdv.size()));
  a(1 + xdv.size()) = 4.0;
  EXPECT_FLOAT_EQ(4.0, a(1 + xdv.size()));
}
