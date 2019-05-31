#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST(test_unit_math_test_ad, to_std_vector) {
  Eigen::VectorXd u(0);
  EXPECT_EQ(0, stan::test::to_std_vector(u).size());

  Eigen::VectorXd x(3);
  x << 1, 2, 3;
  std::vector<double> y = stan::test::to_std_vector(x);
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(x(i), y[i]);
}

TEST(test_unit_math_test_ad, to_eigen_vector) {
  std::vector<double> u;
  EXPECT_EQ(0, stan::test::to_eigen_vector(u).size());

  std::vector<double> v{1, 2, 3};
  Eigen::VectorXd vv = stan::test::to_eigen_vector(v);
  EXPECT_EQ(3, vv.size());
  for (int i = 0; i < 3; ++i)
    EXPECT_EQ(v[i], vv(i));
}

// this is just a dummy to evaluate unary ad tests
template <typename T>
Eigen::Matrix<T, -1, -1> foo(const Eigen::Matrix<T, -1, -1>& x) {
  return -2 * x;
}

TEST(test_unit_math_test_ad, test_ad_unary) {
  Eigen::MatrixXd x(2, 3);
  x << 1, 2, 3, 4, 5, 6;

  // stan::test::expect_ad(foo, x);  // ILLEGAL: can't infer type of foo

  auto g = [](const auto& u) { return foo(u); };
  stan::test::expect_ad(g, x);
}

// dummy to evaluate binary ad tests
// template <typename T1, typename T2>
// Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type, -1,
// -1> bar(const Eigen::Matrix<T1, -1, -1>& x, const T2& y) {
//   using Eigen::Matrix;
//   using boost::math::tools::promote_args;
//   typedef Matrix<typename promote_args<T1, T2>::type, -1, -1> return_t;
//   Eigen::Matrix<return_t, -1, -1> z(x.rows(), x.cols());
//   // for (int i = 0; i < z.size(); ++i) {
//   // return_t a = x(i);
//   // a += y;
//   // z(i) = a;
//   // }
//   return z;
// }

TEST(test_unit_math_test_ad, test_ad_binary) {
  Eigen::MatrixXd x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  double y = -3;

  auto g
      = [](const auto& u, const auto& v) { return stan::math::multiply(u, v); };
  stan::test::expect_ad(g, x, y);
}
