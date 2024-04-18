#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <string>

template <typename T_prob, typename T_prior_sample_size>
void expect_propto_dirichlet_lpdf(T_prob theta, T_prior_sample_size alpha,
                                  T_prob theta2, T_prior_sample_size alpha2,
                                  std::string message) {
  expect_eq_diffs(stan::math::dirichlet_lpdf<false>(theta, alpha),
                  stan::math::dirichlet_lpdf<false>(theta2, alpha2),
                  stan::math::dirichlet_lpdf<true>(theta, alpha),
                  stan::math::dirichlet_lpdf<true>(theta2, alpha2), message);
}

class AgradDistributionsDirichlet : public ::testing::Test {
 protected:
  virtual void SetUp() {
    theta.resize(3, 1);
    theta << 0.2, 0.3, 0.5;
    theta2.resize(3, 1);
    theta2 << 0.1, 0.2, 0.7;
    alpha.resize(3, 1);
    alpha << 1.0, 1.4, 3.2;
    alpha2.resize(3, 1);
    alpha2 << 13.1, 12.9, 10.1;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta;
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha;
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha2;
};

TEST_F(AgradDistributionsDirichlet, Propto) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::to_var;
  using stan::math::var;
  expect_propto_dirichlet_lpdf(to_var(theta), to_var(alpha), to_var(theta2),
                               to_var(alpha2), "var: theta and alpha");
}
TEST_F(AgradDistributionsDirichlet, ProptoTheta) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::to_var;
  using stan::math::var;
  expect_propto_dirichlet_lpdf(to_var(theta), alpha, to_var(theta2), alpha,
                               "var: theta");
}
TEST_F(AgradDistributionsDirichlet, ProptoAlpha) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::to_var;
  using stan::math::var;
  expect_propto_dirichlet_lpdf(theta, to_var(alpha), theta, to_var(alpha2),
                               "var: alpha");
}
TEST_F(AgradDistributionsDirichlet, Bounds) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::dirichlet_lpdf;
  using stan::math::to_var;
  using stan::math::var;
  Matrix<double, Dynamic, 1> good_alpha(2, 1), bad_alpha(2, 1);
  Matrix<double, Dynamic, 1> good_theta(2, 1), bad_theta(2, 1);

  good_theta << 0.25, 0.75;
  good_alpha << 2, 3;
  EXPECT_NO_THROW(dirichlet_lpdf(to_var(good_theta), to_var(good_alpha)));
  EXPECT_NO_THROW(dirichlet_lpdf(to_var(good_theta), good_alpha));
  EXPECT_NO_THROW(dirichlet_lpdf(good_theta, to_var(good_alpha)));

  good_theta << 1.0, 0.0;
  good_alpha << 2, 3;
  EXPECT_NO_THROW(dirichlet_lpdf(to_var(good_theta), to_var(good_alpha)))
      << "elements of theta can be 0";
  EXPECT_NO_THROW(dirichlet_lpdf(to_var(good_theta), good_alpha))
      << "elements of theta can be 0";
  EXPECT_NO_THROW(dirichlet_lpdf(good_theta, to_var(good_alpha)))
      << "elements of theta can be 0";

  bad_theta << 0.25, 0.25;
  EXPECT_THROW(dirichlet_lpdf(to_var(bad_theta), to_var(good_alpha)),
               std::domain_error)
      << "sum of theta is not 1";
  EXPECT_THROW(dirichlet_lpdf(to_var(bad_theta), good_alpha), std::domain_error)
      << "sum of theta is not 1";
  EXPECT_THROW(dirichlet_lpdf(bad_theta, to_var(good_alpha)), std::domain_error)
      << "sum of theta is not 1";

  bad_theta << -0.25, 1.25;
  EXPECT_THROW(dirichlet_lpdf(to_var(bad_theta), to_var(good_alpha)),
               std::domain_error)
      << "theta has element less than 0";
  EXPECT_THROW(dirichlet_lpdf(to_var(bad_theta), good_alpha), std::domain_error)
      << "theta has element less than 0";
  EXPECT_THROW(dirichlet_lpdf(bad_theta, to_var(good_alpha)), std::domain_error)
      << "theta has element less than 0";

  bad_theta << -0.25, 1.25;
  EXPECT_THROW(dirichlet_lpdf(to_var(bad_theta), to_var(good_alpha)),
               std::domain_error)
      << "theta has element less than 0";
  EXPECT_THROW(dirichlet_lpdf(to_var(bad_theta), good_alpha), std::domain_error)
      << "theta has element less than 0";
  EXPECT_THROW(dirichlet_lpdf(bad_theta, to_var(good_alpha)), std::domain_error)
      << "theta has element less than 0";

  bad_alpha << 0.0, 1.0;
  EXPECT_THROW(dirichlet_lpdf(to_var(good_theta), to_var(bad_alpha)),
               std::domain_error)
      << "alpha has element equal to 0";
  EXPECT_THROW(dirichlet_lpdf(to_var(good_theta), bad_alpha), std::domain_error)
      << "alpha has element equal to 0";
  EXPECT_THROW(dirichlet_lpdf(good_theta, to_var(bad_alpha)), std::domain_error)
      << "alpha has element equal to 0";

  bad_alpha << -0.5, 1.0;
  EXPECT_THROW(dirichlet_lpdf(to_var(good_theta), to_var(bad_alpha)),
               std::domain_error)
      << "alpha has element less than 0";
  EXPECT_THROW(dirichlet_lpdf(to_var(good_theta), bad_alpha), std::domain_error)
      << "alpha has element less than 0";
  EXPECT_THROW(dirichlet_lpdf(good_theta, to_var(bad_alpha)), std::domain_error)
      << "alpha has element less than 0";

  bad_alpha = Matrix<double, Dynamic, 1>(4, 1);
  bad_alpha << 1, 2, 3, 4;
  EXPECT_THROW(dirichlet_lpdf(to_var(good_theta), to_var(bad_alpha)),
               std::invalid_argument)
      << "size mismatch: theta is a 2-vector, alpha is a 4-vector";
  EXPECT_THROW(dirichlet_lpdf(to_var(good_theta), bad_alpha),
               std::invalid_argument)
      << "size mismatch: theta is a 2-vector, alpha is a 4-vector";
  EXPECT_THROW(dirichlet_lpdf(good_theta, to_var(bad_alpha)),
               std::invalid_argument)
      << "size mismatch: theta is a 2-vector, alpha is a 4-vector";
}
