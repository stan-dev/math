#include <test/unit/math/test_ad.hpp>
#include <vector>

/**
 supported signatures in stanc3
 gp_dot_prod_cov(array[] real, real) => matrix
 gp_dot_prod_cov(array[] real, array[] real, real) => matrix
 gp_dot_prod_cov(array[] vector, real) => matrix
 gp_dot_prod_cov(array[] vector, array[] vector, real) => matrix
*/
TEST(mathMixMatFun, gp_dot_prod_cov) {
  auto f = [](const auto& x, const auto& sigma) {
    return stan::math::gp_dot_prod_cov(x, sigma);
  };

  auto ff = [](const auto& x1, const auto& x2, const auto& sigma) {
    return stan::math::gp_dot_prod_cov(x1, x2, sigma);
  };

  double sigma = 0.2;
  std::vector<double> x_stdreal{-2, -1, -0.5};
// gp_dot_prod_cov(array[] real, real) => matrix
  stan::test::expect_ad(f, x_stdreal, sigma);
// gp_dot_prod_cov(array[] real, array[] real, real) => matrix
  stan::test::expect_ad(ff, x_stdreal, x_stdreal, sigma);
  Eigen::VectorXd x_vec = Eigen::VectorXd::Random(2);
  std::vector<Eigen::VectorXd> x_stdvec_vec{x_vec, x_vec};
// gp_dot_prod_cov(array[] vector, real) => matrix
  stan::test::expect_ad(f, x_stdvec_vec, sigma);
// gp_dot_prod_cov(array[] vector, array[] vector, real) => matrix
  stan::test::expect_ad(ff, x_stdvec_vec, x_stdvec_vec, sigma);
}
