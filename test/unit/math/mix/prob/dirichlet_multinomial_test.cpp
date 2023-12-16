#include <test/unit/math/test_ad.hpp>

TEST(ProbDistributionsDirichletMultinomial, dirichlet_multinomial) {
  // bind integer vector arg because can't autodiff through
  auto f = [](const std::vector<int>& y) {
    return [=](const auto& alpha) {
      auto lp = stan::math::dirichlet_multinomial_lpmf(y, alpha);
      return lp;
    };
  };

  std::vector<int> y1 = {1, 2, 3, 4};
  std::vector<int> y2 = {30, 75, 409, 34};
  // test if zero counts are handled correctly
  std::vector<int> y3 = {0, 5, 0, 10};

  Eigen::VectorXd alpha1(4);
  alpha1 << 1.0, 2.0, 3.0, 4.0;

  Eigen::VectorXd alpha2(4);
  alpha2 << 1.0, 0.5, 0.1, 1.4;

  // make sure that test cases don't throw exceptions, because expect_ad accepts
  // exceptions, as long as they are thrown consistently

  EXPECT_NO_THROW(f(y1)(alpha1));
  EXPECT_NO_THROW(f(y1)(alpha2));
  EXPECT_NO_THROW(f(y2)(alpha1));
  EXPECT_NO_THROW(f(y2)(alpha2));
  EXPECT_NO_THROW(f(y3)(alpha1));

  stan::test::expect_ad(f(y1), alpha1);
  stan::test::expect_ad(f(y1), alpha2);
  stan::test::expect_ad(f(y2), alpha1);
  stan::test::expect_ad(f(y2), alpha2);
  stan::test::expect_ad(f(y3), alpha1);
}
