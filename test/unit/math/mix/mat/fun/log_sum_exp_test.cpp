#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(MathMixMatFun, logSumExp) {
  auto f = [](const auto& x) { return stan::math::log_sum_exp(x); };

  Eigen::VectorXd x0(0);
  stan::test::expect_ad(f, x0);

  Eigen::VectorXd x1(1);
  x1 << 0;

  Eigen::VectorXd x2(2);
  x2 << 5, 2;

  Eigen::VectorXd x2b(2);
  x2b << 4.9, -std::numeric_limits<double>::infinity();

  Eigen::VectorXd x2c(2);
  x2c << std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity();

  Eigen::VectorXd x4(4);
  x4 << 1, 2, 3, 4;

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  for (const auto& x : std::vector<Eigen::VectorXd>{x1, x2, x2b, x2c, x4}) {
    stan::test::expect_ad(tols, f, x);
    Eigen::RowVectorXd rx = x;
    stan::test::expect_ad(tols, f, rx);
    std::vector<double> stx = std::vector<double>(x.data(), x.data() + x.size());
    stan::test::expect_ad(tols, f, stx);
  }


  Eigen::MatrixXd x23(2, 2);
  x23 << 1, 2, 3, 4;
  stan::test::expect_ad(f, x23);

  std::vector<Eigen::VectorXd> stvx{x2, x2b, x2c};
  stan::test::expect_ad(tols, f, stvx);
  std::vector<Eigen::RowVectorXd> strx{x2.transpose(),
                                       x2b.transpose(),
                                       x2c.transpose()};
  stan::test::expect_ad(tols, f, strx);
  std::vector<std::vector<double>> ststx{std::vector<double>(x2.data(), x2.data() + x2.size()),
                                         std::vector<double>(x2b.data(), x2b.data() + x2b.size()),
                                         std::vector<double>(x2c.data(), x2c.data() + x2c.size())};
  stan::test::expect_ad(tols, f, ststx);
}
