#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_scalar_arguments) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;

  EXPECT_NEAR(-2.12026, poisson_binomial_lcdf(0, p), 0.001);
  EXPECT_NEAR(-0.597837, poisson_binomial_lcdf(1, p), 0.001);
  EXPECT_NEAR(-0.0725707, poisson_binomial_lcdf(2, p), 0.001);
  EXPECT_NEAR(0, poisson_binomial_lcdf(3, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_vectorial_y) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 2};

  EXPECT_NEAR(-0.0725707 * 2, poisson_binomial_lcdf(y, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_vectorial_y_and_theta) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 1};
  std::vector<vec> ps{p, p};

  EXPECT_NEAR(-0.0725707 - 0.597837, poisson_binomial_lcdf(y, ps), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_check_error_scalar_y_oob) {
  static double inff = std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3);
  theta << 0.5, 0.2, 0.7;

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(-1, theta), std::domain_error);
  EXPECT_THROW(stan::math::poisson_binomial_lcdf(4, theta), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_check_error_theta_is_not_prob) {
  static double inff = std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3);
  theta << 0.5, 0.2, 1.1;
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta2(3);
  theta2 << 0.5, 0.2, inff;

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(1, theta), std::domain_error);
  EXPECT_THROW(stan::math::poisson_binomial_lcdf(1, theta2), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_check_error_vectorial_y_oob) {
  static double inff = std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3);
  theta << 0.5, 0.2, 0.1;
  std::vector<int> ys1{-1, 2};
  std::vector<int> ys2{4, 3};

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(ys1, theta),
               std::domain_error);
  EXPECT_THROW(stan::math::poisson_binomial_lcdf(ys2, theta),
               std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     lcdf_check_error_vectorial_y_theta_is_not_prob) {
  static double inff = std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, 1> theta(3);
  theta << inff, 0.2, 0.1;
  std::vector<int> y{0, 2};

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(y, theta), std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     lcdf_check_error_vectorial_theta_is_not_prob) {
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  static double inff = std::numeric_limits<double>::infinity();

  vec theta1(3);
  theta1 << -0.1, 0.2, 0.1;
  vec theta2(3);
  theta2 << 0.5, 0.2, 1.1;
  vec theta3(3);
  theta3 << 0.5, 0.2, inff;

  std::vector<int> ys{0, 2, 3};
  std::vector<vec> thetas{theta1, theta2, theta3};

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(ys, thetas),
               std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     lcdf_check_error_vectorial_y_oob_with_vectorial_theta) {
  using mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec theta1(3);
  theta1 << 0.5, 0.2, 0.1;
  vec theta2(3);
  theta2 << 0.5, 0.2, 0.1;

  std::vector<int> ys{-1, 4};
  std::vector<vec> thetas{theta1, theta2};

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(ys, thetas),
               std::domain_error);
}

TEST(ProbDistributionsPoissonBinomial,
     lcdf_check_error_vectorial_y_and_theta_inconsistent_sizes) {
  using mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec theta(3);
  theta << 0.5, 0.2, 0.1;

  std::vector<int> ys{0, 1};
  std::vector<vec> thetas{theta, theta, theta};

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(ys, thetas),
               std::invalid_argument);
}
