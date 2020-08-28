#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>

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
}

TEST(MathMixMatFun, VarValueSoftmax) {
  using stan::math::var;
  using stan::math::var_value;
  using stan::math::softmax;
  using stan::math::sum;
  Eigen::VectorXd A_val(10);
  for (int i = 0; i < 10; ++i) {
    A_val(i) = static_cast<double>(i);
  }
  // vm is var<matrix> and mv is matrix<var>
  Eigen::Matrix<var, -1, 1> A_mv = A_val;
  Eigen::Matrix<var, -1, 1> A_mv_soft = softmax(A_mv);
  var b_mv = sum(A_mv_soft);
  var final_var = b_mv;
  stan::math::grad(final_var.vi_);
  puts("Vals:");
  std::cout << "\nmat<var>: \n" << A_mv.val();
  std::cout << "\nmat<var> soft: \n" << A_mv_soft.val();
  std::cout << "\nmat<var> sum: \n" << b_mv.val();
  puts("\n\nAdjs:");
  std::cout << "\nmat<var>: \n" << A_mv.adj();
  std::cout << "\nmat<var> soft: \n" << A_mv_soft.adj();
  std::cout << "\nmat<var> sum: \n" << b_mv.adj();
  puts("");
/*  EXPECT_MATRIX_FLOAT_EQ(A_vm.val(), A_mv.val());
  EXPECT_MATRIX_FLOAT_EQ(A_vm.adj(), A_mv.adj());
*/}
