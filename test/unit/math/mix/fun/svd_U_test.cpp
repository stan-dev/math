#include <test/unit/math/test_ad.hpp>
//#include <stan/math/rev/core.hpp>
TEST(MathMixMatFun, svd_U) {
  auto f = [](const auto& x) { return stan::math::svd_U(x); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  Eigen::MatrixXd m11(1, 1);
  m11 << 1.1;
  stan::test::expect_ad(tols, f, m11);

  Eigen::MatrixXd m22(2, 2);
  m22 << 3, -5, 7, 11;
  stan::test::expect_ad(tols, f, m22);
  
  Eigen::MatrixXd m23(2, 3);
  m23 << 3, 5, -7, -11, 13, -17;
  Eigen::MatrixXd m32 = m23.transpose();
  /*stan::math::matrix_v test_arr = m32;
  std::cout << test_arr.val_op() << std::endl;
  std::cout << "----------" << std::endl;
  std::cout << f(test_arr) << std:: endl;*/
  stan::test::expect_ad(tols, f, m23);
  //stan::test::expect_ad(tols, f, m32);

  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 4;
  stan::test::expect_ad(tols, f, a22);
}
