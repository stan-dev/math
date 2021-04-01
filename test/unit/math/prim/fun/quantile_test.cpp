#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

// template <typename T>
// void test_quantile_double() {
//   using stan::math::index_type_t;
//   using stan::math::quantile;

//   T c(1);
//   c[0] = 1.7;
// EXPECT_EQ(quantile(c, 0), 1.7);
// EXPECT_EQ(quantile(c, 1), 1.7);
// EXPECT_EQ(quantile(c, 0.33), 1.7);
// EXPECT_EQ(quantile(c, 0.68), 1.7);

// T v(5);
// v[0] = 1.0;
// v[1] = 2.0;
// v[2] = 3.0;
// v[3] = -1.0;
// v[4] = 19.3;
// EXPECT_EQ(quantile(v, 0), -1.);
// EXPECT_EQ(quantile(v, 1), 19.3);
// EXPECT_EQ(quantile(v, 0.2), 1.);
// std::cout << quantile(v, 0) << std::endl;
// std::cout << quantile(v, 0.1) << std::endl;
// std::cout << quantile(v, 0.8) << std::endl;
// std::cout << quantile(v, 1) << std::endl;
// }

// TEST(MathFunctions, quantileStdVecDouble) {
//   test_quantile_double<std::vector<double> >();
// }

TEST(MathFunctions, bla) {
  using stan::math::quantile;
  std::vector<double> v(5);
  // Eigen::VectorXd v(5);
  // Eigen::RowVectorXd v(5);
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  v[3] = -1.0;
  v[4] = 19.3;

  std::cout << quantile(v, 0.01) << std::endl;
  std::cout << quantile(v, 0.1) << std::endl;
  std::cout << quantile(v, 0.2) << std::endl;
  std::cout << quantile(v, 0.3) << std::endl;
  std::cout << quantile(v, 0.4) << std::endl;
  std::cout << quantile(v, 0.5) << std::endl;
  std::cout << quantile(v, 0.6) << std::endl;
  std::cout << quantile(v, 0.7) << std::endl;
  std::cout << quantile(v, 0.8) << std::endl;
  std::cout << quantile(v, 0.9) << std::endl;
  std::cout << quantile(v, 0.99) << std::endl;
}