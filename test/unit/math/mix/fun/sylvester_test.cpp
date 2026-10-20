#include <test/unit/math/test_ad.hpp>

namespace sylvester_solver_test {

// Real scalar depending on both Re X and Im X, so both carry adjoints
template <typename MatX>
stan::math::var weighted_sum(const MatX& X, const Eigen::MatrixXd& W_re,
                             const Eigen::MatrixXd& W_im) {
  stan::math::var result = 0;
  for (Eigen::Index c = 0; c < X.cols(); ++c) {
    for (Eigen::Index r = 0; r < X.rows(); ++r) {
      result += W_re(r, c) * X(r, c).real() + W_im(r, c) * X(r, c).imag();
    }
  }
  return result;
}

}  // namespace sylvester_solver_test

TEST(mathMixMatFun, sylvester_solver) {
  auto f = [](const auto& a, const auto& b, const auto& c) {
    return stan::math::sylvester_solver(a, b, c);
  };

  // A has a complex eigenvalue pair, so the complex Schur path is used
  Eigen::MatrixXd A(2, 2);
  A << 1, -1, 1, 1;
  Eigen::MatrixXd B(3, 3);
  B << 5, 0, 1, 2, 6, 0, -1, 1, 4;
  Eigen::MatrixXd C(2, 3);
  C << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, A, B, C);
  // expect_ad_matvar does not apply: the result has complex<var> scalars,
  // which have no var_value<Matrix> form. See sylvester_solver_varmat.

  // Lyapunov: B = A^T shares A's autodiff variables
  auto g = [](const auto& a, const auto& c) {
    return stan::math::sylvester_solver(a, a.transpose(), c);
  };
  Eigen::MatrixXd Q(2, 2);
  Q << 1, 2, 3, 4;
  stan::test::expect_ad(g, A, Q);

  // size mismatch must throw for every autodiff type, as it does for double
  Eigen::MatrixXd C_bad(3, 3);
  C_bad.setOnes();
  EXPECT_THROW(stan::math::sylvester_solver(A, B, C_bad),
               std::invalid_argument);

  auto A_v = stan::math::to_var(A);
  auto B_v = stan::math::to_var(B);
  auto C_bad_v = stan::math::to_var(C_bad);
  EXPECT_THROW(stan::math::sylvester_solver(A_v, B_v, C_bad_v),
               std::invalid_argument);
  stan::math::recover_memory();
}

TEST(mathMixMatFun, sylvester_solver_complex) {
  auto f = [](const auto& a, const auto& b, const auto& c) {
    return stan::math::sylvester_solver(a, b, c);
  };
  using cdouble = std::complex<double>;

  // eigenvalue sums lambda_i(A) + mu_j(B) have real part around 9
  Eigen::MatrixXcd A(2, 2);
  A << cdouble(4, 1), cdouble(1, -0.5), cdouble(-1, 0.3), cdouble(3, -0.2);
  Eigen::MatrixXcd B(2, 2);
  B << cdouble(5, -0.4), cdouble(0.5, 1), cdouble(2, -0.7), cdouble(6, 0.3);
  Eigen::MatrixXcd C(2, 2);
  C << cdouble(1, 2), cdouble(-0.5, 1), cdouble(3, -1), cdouble(0.7, 0.2);
  stan::test::expect_ad(f, A, B, C);

  // real and complex arguments mixed
  Eigen::MatrixXd A_r = A.real();
  Eigen::MatrixXd B_r = B.real();
  // real system, complex right-hand side
  stan::test::expect_ad(f, A_r, B_r, C);
}

TEST(mathMixMatFun, sylvester_solver_varmat) {
  // var_value<Matrix> inputs must give the same value and adjoints as
  // Matrix<var> inputs
  using stan::math::var;
  using stan::math::var_value;
  using sylvester_solver_test::weighted_sum;
  Eigen::MatrixXd A(2, 2);
  A << 4, 1, -1, 3;
  Eigen::MatrixXd B(3, 3);
  B << 5, 0, 1, 2, 6, 0, -1, 1, 4;
  Eigen::MatrixXd C(2, 3);
  C << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd W_re(2, 3);
  W_re << 0.3, -1.2, 0.7, 2.1, -0.4, 1.5;
  Eigen::MatrixXd W_im(2, 3);
  W_im << -0.8, 0.6, 1.9, -1.1, 0.2, -0.5;

  var_value<Eigen::MatrixXd> A_vm(A);
  var_value<Eigen::MatrixXd> B_vm(B);
  var_value<Eigen::MatrixXd> C_vm(C);
  var f_vm = weighted_sum(stan::math::sylvester_solver(A_vm, B_vm, C_vm), W_re,
                          W_im);
  f_vm.grad();
  const double val_vm = f_vm.val();
  Eigen::MatrixXd adj_A_vm = A_vm.adj();
  Eigen::MatrixXd adj_B_vm = B_vm.adj();
  Eigen::MatrixXd adj_C_vm = C_vm.adj();
  stan::math::recover_memory();

  auto A_mv = stan::math::to_var(A);
  auto B_mv = stan::math::to_var(B);
  auto C_mv = stan::math::to_var(C);
  var f_mv = weighted_sum(stan::math::sylvester_solver(A_mv, B_mv, C_mv), W_re,
                          W_im);
  f_mv.grad();

  stan::test::expect_near_rel("value", val_vm, f_mv.val());
  stan::test::expect_near_rel("adjoint of A", adj_A_vm,
                              Eigen::MatrixXd(A_mv.adj()));
  stan::test::expect_near_rel("adjoint of B", adj_B_vm,
                              Eigen::MatrixXd(B_mv.adj()));
  stan::test::expect_near_rel("adjoint of C", adj_C_vm,
                              Eigen::MatrixXd(C_mv.adj()));
  stan::math::recover_memory();
}

TEST(mathMixMatFun, sylvester_solver_result_is_independent) {
  // Writing into the result must not change what the reverse pass reads:
  // after X(0, 0) = 2 X(0, 0), d Re X(0, 0) / dC must double
  using stan::math::var;
  Eigen::MatrixXd A(2, 2);
  A << 4, 1, -1, 3;
  Eigen::MatrixXd B(2, 2);
  B << 5, 0, 2, 6;
  Eigen::MatrixXd C(2, 2);
  C << 1, 2, 3, 4;

  auto C1 = stan::math::to_var(C);
  auto X1 = stan::math::sylvester_solver(A, B, C1);
  var f1 = X1(0, 0).real();
  f1.grad();
  Eigen::MatrixXd adj_C = C1.adj();
  stan::math::recover_memory();

  auto C2 = stan::math::to_var(C);
  auto X2 = stan::math::sylvester_solver(A, B, C2);
  X2(0, 0) = 2.0 * X2(0, 0);
  var f2 = X2(0, 0).real();
  f2.grad();
  stan::test::expect_near_rel("adjoint of C", Eigen::MatrixXd(2.0 * adj_C),
                              Eigen::MatrixXd(C2.adj()));
  stan::math::recover_memory();
}
