#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_ordered_logistic_glm_lpmf) {
  auto f = [](const auto y) {
    return [=](const auto& x, const auto& beta, const auto& cutpoints) {
      return stan::math::ordered_logistic_glm_lpmf(y, x, beta, cutpoints);
    };
  };

  std::vector<int> y{0, 1, 2};
  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);
  Eigen::RowVectorXd x_rowvec = x.row(0);
  Eigen::VectorXd cutpoints = stan::math::sort_asc(Eigen::VectorXd::Random(2));
  Eigen::VectorXd beta = Eigen::VectorXd::Random(2);

  stan::test::expect_ad(f(y[0]), x, beta, cutpoints);
  stan::test::expect_ad(f(y[0]), x_rowvec, beta, cutpoints);
  stan::test::expect_ad(f(y), x, beta, cutpoints);
  stan::test::expect_ad(f(y), x_rowvec, beta, cutpoints);
}
