#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, solve_lyapunov) {
  auto f = [](const auto& a, const auto& c) {
    return stan::math::solve_lyapunov(a, c);
  };

  // A has a complex eigenvalue pair, so the complex Schur path is used
  Eigen::MatrixXd A(3, 3);
  A << -1, -1, 0, 1, -1, 0.5, 0.3, -0.2, -2;
  Eigen::MatrixXd C(3, 3);
  C << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::test::expect_ad(f, A, C);

  // size mismatch must throw for every autodiff type, as it does for double
  Eigen::MatrixXd C_bad(3, 2);
  C_bad.setOnes();
  EXPECT_THROW(stan::math::solve_lyapunov(A, C_bad), std::invalid_argument);

  auto A_v = stan::math::to_var(A);
  auto C_bad_v = stan::math::to_var(C_bad);
  EXPECT_THROW(stan::math::solve_lyapunov(A_v, C_bad_v), std::invalid_argument);
  stan::math::recover_memory();
}

TEST(mathMixMatFun, solve_lyapunov_complex) {
  auto f = [](const auto& a, const auto& c) {
    return stan::math::solve_lyapunov(a, c);
  };
  using cdouble = std::complex<double>;

  // eigenvalue sums lambda_i + conj(lambda_j) have real part around -7
  Eigen::MatrixXcd A(2, 2);
  A << cdouble(-4, 1), cdouble(1, -0.5), cdouble(-1, 0.3), cdouble(-3, -0.2);
  Eigen::MatrixXcd C(2, 2);
  C << cdouble(1, 2), cdouble(-0.5, 1), cdouble(3, -1), cdouble(0.7, 0.2);
  stan::test::expect_ad(f, A, C);

  // real system, complex right-hand side
  Eigen::MatrixXd A_r = A.real();
  stan::test::expect_ad(f, A_r, C);
}

TEST(mathMixMatFun, solve_lyapunov_triangular) {
  // Schur of a triangular matrix deflates on exact zeros, which must not
  // drop the tangents of the zero entries
  auto f = [](const auto& a, const auto& c) {
    return stan::math::solve_lyapunov(a, c);
  };
  Eigen::MatrixXd A(3, 3);
  A << -5, 2, -1, 0, -6, 1, 0, 0, -4;
  Eigen::MatrixXd C(3, 3);
  C << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::test::expect_ad(f, A, C);
}
