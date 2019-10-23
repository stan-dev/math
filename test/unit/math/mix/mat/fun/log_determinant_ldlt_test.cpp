#include <test/unit/math/test_ad.hpp>

template <typename T, int R, int C>
stan::math::LDLT_factor<T, R, C> ldlt_factor(const Eigen::Matrix<T, R, C>& x) {
  stan::math::LDLT_factor<T, R, C> ldlt;
  ldlt.compute(x);
  return ldlt;
}

TEST(MathMixMatFun, logDeterminantLdlt) {
  auto f = [](const auto& x) {
    auto z = ((x + x.transpose()) * 0.5).eval();
    auto y = ldlt_factor(z);
    return stan::math::log_determinant_ldlt(y);
  };

  Eigen::MatrixXd a(2, 2);
  a << 3, 0, 0, 4;
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd c(2, 2);
  c << 1, 0, 0, 3;
  stan::test::expect_ad(f, c);

  Eigen::MatrixXd b(2, 2);
  b << 2, 1, 1, 3;
  stan::test::expect_ad(f, b);

  tols.gradient_grad_ = 5;
  for (const auto& rho : std::vector<double>{0, 0.9}) {
    for (const auto& y : stan::test::ar_test_cov_matrices(1, 3, rho)) {
      stan::test::expect_ad(f, y);
    }
  }
}
