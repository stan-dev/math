#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, softmax) {
  auto f = [](const auto& x) { return stan::math::softmax(x); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  Eigen::VectorXd a(0);
  stan::test::expect_ad(tols, f, a);

  Eigen::VectorXd b(1);
  b << 0;
  stan::test::expect_ad(tols, f, b);

  Eigen::VectorXd c(2);
  c << -1, 1;
  stan::test::expect_ad(tols, f, c);

  Eigen::VectorXd d(3);
  d << -1, 1, 10;
  stan::test::expect_ad(tols, f, d);

  Eigen::VectorXd d2(3);
  d2 << 0.5, -1, 3;
  stan::test::expect_ad(tols, f, d2);

  Eigen::VectorXd d3(3);
  d3 << 4, 3, -2;
  stan::test::expect_ad(tols, f, d3);

  Eigen::VectorXd d4(3);
  d4 << 0, 3, -1;
  stan::test::expect_ad(tols, f, d4);

  // Row Vectors
  Eigen::RowVectorXd rx0(0);  // error case
  stan::test::expect_ad(tols, f, rx0);

  Eigen::RowVectorXd rx1(1);
  rx1 << 0;
  stan::test::expect_ad(tols, f, rx1);

  Eigen::RowVectorXd rx2(2);
  rx2 << -1, 1;
  stan::test::expect_ad(tols, f, rx2);

  Eigen::RowVectorXd rx3(3);
  rx3 << -1, 1, 10;
  stan::test::expect_ad(tols, f, rx3);

  Eigen::RowVectorXd rx3b(3);
  rx3b << 0, 1, 2;
  stan::test::expect_ad(tols, f, rx3b);

  Eigen::RowVectorXd rx3c(3);
  rx3c << 2, 1, 1;
  stan::test::expect_ad(tols, f, rx3c);

  // std vectors
  std::vector<double> stx0(0);  // error case
  stan::test::expect_ad(tols, f, stx0);

  std::vector<double> stx1{0};
  stan::test::expect_ad(tols, f, stx1);

  std::vector<double> stx2{-1, 1};
  stan::test::expect_ad(tols, f, stx2);

  std::vector<double> stx3{-1, 1, 10};
  stan::test::expect_ad(tols, f, stx3);

  std::vector<double> stx3b{0, 1, 2};
  stan::test::expect_ad(tols, f, stx3b);

  std::vector<double> stx3c{2, 1, 1};
  stan::test::expect_ad(tols, f, stx3c);

  // Nested containers
  std::vector<Eigen::VectorXd> stvx0{a, a};  // error case
  stan::test::expect_ad(tols, f, stvx0);

  std::vector<Eigen::VectorXd> stvx1{b, b};
  stan::test::expect_ad(tols, f, stvx1);

  std::vector<Eigen::RowVectorXd> strx0{rx0, rx0};  // error case
  stan::test::expect_ad(tols, f, strx0);

  std::vector<Eigen::RowVectorXd> strx1{rx1, rx1};
  stan::test::expect_ad(tols, f, strx1);

  std::vector<std::vector<double>> ststx0{stx0, stx0};  // error case
  stan::test::expect_ad(tols, f, ststx0);

  std::vector<std::vector<double>> ststx1{stx1, stx1};
  stan::test::expect_ad(tols, f, ststx1);
}
