#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(MathMixMatFun, qr) {
  auto f = [](const auto& x) { return std::get<0>(stan::math::qr(x)); };
  auto g = [](const auto& x) { return std::get<1>(stan::math::qr(x)); };

  Eigen::MatrixXd a(0, 0);
  stan::test::expect_ad(f, a);
  stan::test::expect_ad(g, a);

  Eigen::MatrixXd b(1, 1);
  b << 1.5;
  stan::test::expect_ad(f, b);
  stan::test::expect_ad(g, b);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  Eigen::MatrixXd c(3, 2);
  c << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(tols, f, c);
  stan::test::expect_ad(tols, g, c);

  // cols > rows case
  Eigen::MatrixXd b_tr = b.transpose();
  stan::test::expect_ad(f, b_tr);
  stan::test::expect_ad(g, b_tr);
}
