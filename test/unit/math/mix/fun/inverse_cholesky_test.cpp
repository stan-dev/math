#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include <limits>

namespace inverse_cholesky_test {
// can't autodiff directly through Cholesky due to symmetry test;
// use unconstrained input and constrain to test Cholesky derivs;
// dof must be (n choose 2) + n
using stan::math::cholesky_decompose;

auto f(int dof) {
  return [=](const auto& x) {
    stan::math::promote_scalar_t<stan::value_type_t<decltype(x)>,
                                 Eigen::Matrix<double, -1, -1>>
        y;
    try {
      y = stan::math::cov_matrix_constrain(x, dof);
    } catch (...) {
      ADD_FAILURE() << "FAILED AT COV_MATRIX_CONSTRAIN";
      throw;
    }
    return stan::math::inverse_cholesky(cholesky_decompose(y));
  };
}

auto f_matvar = [](const auto& x) { return stan::math::inverse_cholesky(x); };

}  // namespace cholesky_decompose_test

TEST(MathMixMatFun, inverseCholeskyDecomposeSpecific) {
  using stan::math::cholesky_decompose;
  // 1 x 1 matrix;  (1 choose 2) + 1 = 1
  Eigen::VectorXd x1(1);
  x1 << 1;
  stan::test::expect_ad(inverse_cholesky_test::f(1), x1);
  Eigen::MatrixXd x1_mat = cholesky_decompose(stan::math::cov_matrix_constrain(x1, 1));
 // stan::test::expect_ad_matvar(inverse_cholesky_test::f_matvar, x1_mat);

  Eigen::MatrixXd j(3, 3);
  j << 1, 0, 0, 2, 5, 0, 1, 2, 3;
  stan::test::expect_ad(inverse_cholesky_test::f_matvar, j);
  stan::test::expect_ad_matvar(inverse_cholesky_test::f_matvar, j);


  // // 2 x 2 matrix;  (2 choose 2) + 2 = 3
   Eigen::VectorXd x3(3);
   x3 << 1, 2, -1;
   stan::test::expect_ad(inverse_cholesky_test::f(2), x3);
  // Eigen::MatrixXd x3_mat = stan::math::cov_matrix_constrain(x3, 2);
  // stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, x3_mat);

  // // 3 x 3 matrix;  (3 choose 2) + 3 = 6
  // Eigen::VectorXd x6(6);
  // x6 << 1, -1, 1.1, 1.4, 2.1, 0.7;
  // stan::test::expect_ad(cholesky_decompose_test::f(3), x6);
  // Eigen::MatrixXd x6_mat = stan::math::cov_matrix_constrain(x6, 3);
  // stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, x6_mat);

  // // 4 x 4 matrix;  (4 choose 2) + 4 = 10
  // Eigen::VectorXd x10(10);
  // x10 << 1, -0.1, 1.1, 1.4, -1.1, 0.7, 1.0, 1.3, -0.5, 0.3;
  // stan::test::expect_ad(cholesky_decompose_test::f(4), x10);
  // Eigen::MatrixXd x10_mat = stan::math::cov_matrix_constrain(x10, 4);
  // stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, x10_mat);

  // // 2 x 3 matrix will throw; test directly
  // auto g = [](const auto& x) { return stan::math::cholesky_decompose(x); };
  // Eigen::MatrixXd y(2, 3);
  // y << 1, 2, 3, 4, 5, 6;
  // stan::test::expect_ad(g, y);
  // stan::test::expect_ad_matvar(g, y);

  // // asymmetric will throw
  // Eigen::MatrixXd z(2, 2);
  // z << 1, 2, 3, 4;
  // stan::test::expect_ad(g, z);
  // stan::test::expect_ad_matvar(g, y);
}

// TEST(MathMixMatFun, inverse_cholesky) {
//   using stan::test::expect_ad;
 
//   auto f = [](const auto& x) {
//     return stan::math::inverse_cholesky(x);
//   };

//   Eigen::MatrixXd m00(0, 0);
//   expect_ad(f, m00);

//   Eigen::MatrixXd m11(1, 1);
//   m11 << 1.3;
//   expect_ad(f, m11);

//   Eigen::MatrixXd a(2, 2);
//   a << 1.414214, 0,  2.121320, 1.581139;
//   expect_ad(f, a);

//   Eigen::MatrixXd b(2, 2);
//   b << 1, -1, -1, -1;
//   expect_ad(f, b);

//   Eigen::MatrixXd c(2, 2);
//   c << 2, 3, 1, 7;
//   expect_ad(f, c);

//   Eigen::MatrixXd j(3, 3);
//   j << 1, 0, 0, 2, 5, 0, 1, 2, 3;
//   expect_ad(f, j);
// }

