#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(MathMixMatFun, kronecker) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::kronecker(x, y);
  };

  Eigen::MatrixXd a23{{1, -2, 3}, {0.5, 4, -1}};
  Eigen::MatrixXd b32{{2, 1}, {0, -3}, {5, 0.25}};
  Eigen::MatrixXd c11{{-1.5}};
  Eigen::MatrixXd v31{{0.5}, {-1}, {2}};
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m02(0, 2);
  Eigen::MatrixXd m20(2, 0);
  for (const auto& x : {a23, b32, c11, v31, m00, m02, m20}) {
    for (const auto& y : {a23, b32, c11, v31, m00, m02, m20}) {
      stan::test::expect_ad(f, x, y);
      stan::test::expect_ad_matvar(f, x, y);
    }
  }

  Eigen::VectorXd v{{0.5, -1, 2}};
  Eigen::RowVectorXd rv{{3, -0.25}};
  stan::test::expect_ad(f, v, rv);
  stan::test::expect_ad_matvar(f, v, rv);
  stan::test::expect_ad(f, rv, a23);
  stan::test::expect_ad_matvar(f, rv, a23);
}

TEST(MathMixMatFun, kronecker_complex) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::kronecker(x, y);
  };
  using c = std::complex<double>;
  Eigen::MatrixXcd a{{c(1, 2), c(-1, 0)}, {c(0, 0.5), c(3, 0)}};
  Eigen::MatrixXcd b{{c(2, -1)}, {c(1, 1)}, {c(0, -2)}};
  stan::test::expect_ad(f, a, b);
  stan::test::expect_ad(f, b, a);
}
