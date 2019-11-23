#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, inverseSpd) {
  auto f = [](const auto& x) {
    if (x.rows() != x.cols())
      return stan::math::inverse_spd(x);
    auto y = ((x + x.transpose()) * 0.5).eval();  // symmetry for finite diffs
    return stan::math::inverse_spd(y);
  };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  Eigen::MatrixXd m11(1, 1);
  m11 << 1.3;
  stan::test::expect_ad(f, m11);

  Eigen::MatrixXd a(2, 2);
  a << 2, 3, 3, 7;
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd b(2, 2);
  b << 1, -1, -1, -1;
  stan::test::expect_ad(f, b);

  Eigen::MatrixXd c(2, 2);
  c << 2, 3, 1, 7;
  stan::test::expect_ad(f, c);

  for (int k = 1; k < 4; ++k) {
    Eigen::MatrixXd d(k, k);
    for (int i = 0; i < k; ++i) {
      d(i, i) = 1;
      for (int j = 0; j < i; ++j) {
        d(i, j) = std::pow(0.9, std::fabs(i - j));
        d(j, i) = d(i, j);
      }
    }
    stan::test::expect_ad(f, d);
  }

  Eigen::MatrixXd h(3, 3);  // not positive definite
  h << 1, 2, 3, 2, 4, 5, 3, 5, 6;
  stan::test::expect_ad(f, h);

  Eigen::MatrixXd j(3, 3);
  j << 2, -1, 0, -1, 2, -1, 0, -1, 2;
  stan::test::expect_ad(f, j);

  // test functor that doesn't symmetrize input to test errors
  auto f_asym = [](const auto& x) { return stan::math::inverse_spd(x); };

  Eigen::MatrixXd g(3, 3);  // not symmetric
  g << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::test::expect_ad(f_asym, g);

  Eigen::MatrixXd e(3, 2);  // not square
  e << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f_asym, e);
}
