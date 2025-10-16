#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, squaredDistance) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::squared_distance(x1, x2);
  };
  stan::test::expect_ad(f, 1.0, 1.0);
  stan::test::expect_ad(f, 1.0, 4.0);
  stan::test::expect_ad(f, 4.0, 1.0);
  stan::test::expect_ad(f, 4.0, 4.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}
void expect_squared_distance(const std::vector<double>& sv1,
                             const std::vector<double>& sv2,
                             const stan::test::ad_tolerances& tols
                             = stan::test::ad_tolerances()) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::squared_distance(x, y);
  };

  Eigen::VectorXd v1 = stan::test::to_vector(sv1);
  Eigen::VectorXd v2 = stan::test::to_vector(sv2);
  Eigen::RowVectorXd rv1 = stan::test::to_row_vector(sv1);
  Eigen::RowVectorXd rv2 = stan::test::to_row_vector(sv2);

  stan::test::expect_ad(tols, f, v1, v2);
  stan::test::expect_ad(tols, f, rv1, rv2);
  stan::test::expect_ad(tols, f, v1, rv2);
  stan::test::expect_ad(tols, f, rv1, v2);

  stan::test::expect_ad_matvar(tols, f, v1, v2);
  stan::test::expect_ad_matvar(tols, f, rv1, rv2);
  stan::test::expect_ad_matvar(tols, f, v1, rv2);
  stan::test::expect_ad_matvar(tols, f, rv1, v2);
}

TEST(MathMixMatFun, squaredDistance) {
  std::vector<double> a0;
  std::vector<double> b0;
  expect_squared_distance(a0, b0);

  std::vector<double> a1{1.1};
  std::vector<double> b1{-2.2};
  expect_squared_distance(a1, b1);

  std::vector<double> a2{1.1, 3.9};
  std::vector<double> b2{-2.2, 15.2};
  expect_squared_distance(a2, b2);

  std::vector<double> a3{1, 3, -5};
  std::vector<double> b3{4, -2, -1};
  expect_squared_distance(a3, b3);

  std::vector<double> c3{-1, 0, 1};
  std::vector<double> d3{1, 2, 3};
  expect_squared_distance(c3, d3);

  std::vector<double> g3{1, 2, 3};
  std::vector<double> h3{4, 5, 6};
  expect_squared_distance(g3, h3);

  // infinite input
  std::vector<double> e3{1, 1, std::numeric_limits<double>::infinity()};
  expect_squared_distance(c3, e3);
  expect_squared_distance(e3, c3);

  // exceptional input
  std::vector<double> a3nan{1, std::numeric_limits<double>::quiet_NaN(), 2};
  expect_squared_distance(a3nan, b3);
  expect_squared_distance(b3, a3nan);

  expect_squared_distance(a2, b3);
  expect_squared_distance(a3, b2);
}
