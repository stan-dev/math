#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(ProbDistributionsPoissonBinomial, lcdf_length_0_length_1) {
  Eigen::VectorXd v0(0);
  Eigen::VectorXd v1(1);
  v1 << 0.4;

  EXPECT_FLOAT_EQ(stan::math::poisson_binomial_lcdf(0, v0), 0.0);
  EXPECT_FLOAT_EQ(stan::math::poisson_binomial_lcdf(1, v1), 0.0);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_length_0_length_1_vectorial_y) {
  Eigen::VectorXd v0(0);
  Eigen::VectorXd v1(1);
  v1 << 0.4;

  std::vector<int> y0{0, 0};
  std::vector<int> y1{1, 1};

  EXPECT_FLOAT_EQ(stan::math::poisson_binomial_lcdf(y0, v0), 0.0);
  EXPECT_FLOAT_EQ(stan::math::poisson_binomial_lcdf(y1, v1), 0.0);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_scalar_arguments) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;

  EXPECT_NEAR(-2.12026, poisson_binomial_lcdf(0, p), 1e-5);
  EXPECT_NEAR(-0.597837, poisson_binomial_lcdf(1, p), 1e-6);
  EXPECT_NEAR(-0.0725707, poisson_binomial_lcdf(2, p), 1e-6);
  EXPECT_NEAR(0, poisson_binomial_lcdf(3, p), 1e-6);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_vectorial_y) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 2};

  EXPECT_NEAR(-0.0725707 * 2, poisson_binomial_lcdf(y, p), 1e-6);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_vectorial_y_and_theta) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 1};
  std::vector<vec> ps{p, p};

  EXPECT_NEAR(-0.0725707 - 0.597837, poisson_binomial_lcdf(y, ps), 1e-6);
}

TEST(ProbDistributionsPoissonBinomial,
     lcdf_works_on_scalar_y_and_vectorial_theta) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  int y = 2;
  std::vector<vec> ps{p};

  EXPECT_NEAR(-0.0725707, poisson_binomial_lcdf(y, ps), 1e-6);
}

TEST(ProbDistributionsPoissonBinomial,
     lcdf_works_on_scalar_y_and_vectorial_theta2) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  int y = 2;
  std::vector<vec> ps{p, p};

  EXPECT_NEAR(-0.0725707 * 2.0, poisson_binomial_lcdf(y, ps), 1e-6);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_check_error_scalar_y_oob) {
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
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec theta(3);
  theta << 0.5, 0.2, 0.1;

  std::vector<int> ys{0, 1};
  std::vector<vec> thetas{theta, theta, theta};

  EXPECT_THROW(stan::math::poisson_binomial_lcdf(ys, thetas),
               std::invalid_argument);
}

TEST(ProbDistributionsPoissonBinomial, log_cdf_matches_lcdf) {
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  int y = 2;
  vec theta(3, 1);
  theta << 0.5, 0.2, 0.7;

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lcdf(y, theta)),
                  (stan::math::poisson_binomial_cdf_log(y, theta)));
}

TEST(ProbDistributionsPoissonBinomial, log_cdf_matches_lcdf_vectorial) {
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  int y = 2;
  vec theta(3, 1);
  theta << 0.5, 0.2, 0.7;
  std::vector<vec> thetas{theta, theta};

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lcdf(y, thetas)),
                  (stan::math::poisson_binomial_cdf_log(y, thetas)));
}
