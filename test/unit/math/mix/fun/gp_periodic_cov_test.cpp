#include <test/unit/math/test_ad.hpp>
#include <vector>


TEST(mathMixMatFun, gp_periodic_cov) {
  auto f = [](const auto& x, const auto& sigma, const auto& l) {
    return stan::math::gp_periodic_cov(x, sigma, l, l);
  };

  double sigma = 0.5;
  double l = 1;
  std::vector<double> x_stdreal{-2, -1, -0.5};
  stan::test::expect_ad(f, x_stdreal, sigma, l);
  Eigen::VectorXd x_vec = Eigen::VectorXd::Random(2);
  std::vector<Eigen::VectorXd> x_stdvec_vec{x_vec, x_vec};
  stan::test::expect_ad(f, x_stdvec_vec, sigma, l);
}

TEST(mathMixMatFun, gp_periodic_cov_stdvecvec) {
  auto ff = [](const auto& x1, const auto& sigma, const auto& l) {
    return stan::math::gp_periodic_cov(x1, x1, sigma, l, l);
  };

  double sigma = 0.5;
  double l = 1;
  std::vector<double> x_stdreal{-2, -1, -0.5};
  stan::test::expect_ad(ff, x_stdreal, sigma, l);
  Eigen::VectorXd x_vec = Eigen::VectorXd::Random(2);
  std::vector<Eigen::VectorXd> x_stdvec_vec{x_vec, x_vec};
  //stan::test::expect_ad(ff, x_stdvec_vec, sigma, l);
}
