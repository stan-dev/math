#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, ProbDistributionsGaussCopulaCholesky) {
  // Bind functors for compatibility with AD framework
  auto f = [](const auto func1, const auto func2) {
    return [=](const auto& y, const auto& args, const auto& sigma) {
      auto lcdf_functors = std::make_tuple(std::make_tuple(func1, args[0]),
                                           std::make_tuple(func2, args[1]));
      auto sigma_sym = stan::math::multiply(0.5, sigma + sigma.transpose());
      auto L = stan::math::cholesky_decompose(sigma_sym);
      return stan::math::gaussian_copula_cholesky_lpdf(y, lcdf_functors, L);
    };
  };

  // y[0] ~ gamma(2, 1)
  // y[1] ~ std_normal()
  Eigen::VectorXd y1(2);
  y1 << 1.0, 0.1;
  Eigen::VectorXd args(2);
  args << 0, 2.0;

  auto func1 = [](const auto& y, const auto& mu) {
    return stan::math::normal_lcdf(y, mu, 1.0);
  };
  auto func2 = [](const auto& y) { return stan::math::std_normal_lcdf(y); };

  Eigen::MatrixXd Sigma22(2, 2);
  Sigma22 << 2.0, 0.5, 0.5, 1.1;

  stan::test::expect_ad(f(func2, func2), y1, args, Sigma22);
}
