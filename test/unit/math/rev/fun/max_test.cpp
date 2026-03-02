#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(MathRev, max_scalar_permutations) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::max(x, y);
  };
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 2.0, 1.0);
  stan::test::expect_ad(f, -5.0, 0.0);
}

TEST(MathRev, max_scalar_tie) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::max(x, y);
  };
  stan::test::expect_ad(f, 3.0, 3.0);
}

TEST(MathRev, max_scalar_boundary) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::max(x, y);
  };
  double nan = std::numeric_limits<double>::quiet_NaN();
  double inf = std::numeric_limits<double>::infinity();
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, inf, -inf);
}

TEST(MathRev, max_container_types) {
  auto f = [](const auto& x) {
    return stan::math::max(x);
  };
  Eigen::VectorXd v(3); v << 1, 2, 3;
  Eigen::RowVectorXd rv(3); rv << 1, 2, 3;
  std::vector<double> sv{1, 2, 3};
  Eigen::MatrixXd m(2, 2); m << 1, 5, 2, 3;

  stan::test::expect_ad(f, v);
  stan::test::expect_ad(f, rv);
  stan::test::expect_ad(f, sv);
  stan::test::expect_ad(f, m);
}

TEST(MathRev, max_container_ties) {
  auto f = [](const auto& x) {
    return stan::math::max(x);
  };
  Eigen::VectorXd v(4); v << 1.5, 1.5, 1.5, 1.5;
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-3; 
  stan::test::expect_ad(tols, f, v);
}

TEST(MathRev, max_container_edges) {
  Eigen::VectorXd v_single(1); v_single << 42.0;
  EXPECT_FLOAT_EQ(42.0, stan::math::value_of(stan::math::max(v_single)));
  
  Eigen::VectorXd v_empty(0);
  EXPECT_FLOAT_EQ(stan::math::NEGATIVE_INFTY, stan::math::max(v_empty));
}

TEST(MathRev, max_container_nan) {
  auto f = [](const auto& x) {
    return stan::math::max(x);
  };

  double nan = std::numeric_limits<double>::quiet_NaN();
  
  Eigen::VectorXd v1(3);
  v1 << 1.0, nan, 2.0;
  stan::test::expect_ad(f, v1);

  Eigen::VectorXd v2(2);
  v2 << nan, nan;
  stan::test::expect_ad(f, v2);
}

TEST(mathRevFun, max_fvar_check) {
  using stan::math::fvar;
  using stan::math::max;

  fvar<double> x = 1.0;
  fvar<double> y = 2.0;
  
  EXPECT_NO_THROW(max(x, y));
}