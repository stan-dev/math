#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(MathRev, min_scalar_permutations) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::min(x, y);
  };
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 2.0, 1.0);
}

TEST(MathRev, min_scalar_tie) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::min(x, y);
  };
  stan::test::expect_ad(f, 0.0, 0.0);
}

TEST(MathRev, min_scalar_boundary) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::min(x, y);
  };
  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, 1.0, nan);
}

TEST(MathRev, min_container_types) {
  auto f = [](const auto& x) {
    return stan::math::min(x);
  };
  Eigen::VectorXd v(3); v << 3, 2, 1;
  Eigen::RowVectorXd rv(3); rv << 3, 2, 1;
  std::vector<double> sv{3, 2, 1};
  Eigen::MatrixXd m(2, 2); m << 1, 0.5, 2, 3;

  stan::test::expect_ad(f, v);
  stan::test::expect_ad(f, rv);
  stan::test::expect_ad(f, sv);
  stan::test::expect_ad(f, m);
}

TEST(MathRev, min_container_ties) {
  auto f = [](const auto& x) {
    return stan::math::min(x);
  };
  Eigen::VectorXd v(4); v << -1.1, -1.1, -1.1, -1.1;
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-3; 
  stan::test::expect_ad(tols, f, v);
}

TEST(MathRev, min_container_edges) {
  Eigen::VectorXd v_single(1); v_single << -42.0;
  EXPECT_FLOAT_EQ(-42.0, stan::math::value_of(stan::math::min(v_single)));
  
  Eigen::VectorXd v_empty(0);
  EXPECT_FLOAT_EQ(stan::math::INFTY, stan::math::min(v_empty));
}

TEST(MathRev, min_container_nan) {
  auto f = [](const auto& x) {
    return stan::math::min(x);
  };

  double nan = std::numeric_limits<double>::quiet_NaN();
  
  Eigen::VectorXd v1(3);
  v1 << 10.0, nan, -5.0;
  stan::test::expect_ad(f, v1);

  Eigen::VectorXd v2(2);
  v2 << nan, nan;
  stan::test::expect_ad(f, v2);
}