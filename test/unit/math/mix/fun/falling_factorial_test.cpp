#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fallingFactorial) {
  auto f = [](const int x2) {
    return
        [=](const auto& x1) { return stan::math::falling_factorial(x1, x2); };
  };
  stan::test::expect_ad(f(-2), -3.0);  // throws

  stan::test::expect_ad(f(3), 5);

  // 3rd order derivatives near zero
  stan::test::ad_tolerances tols;
  tols.grad_hessian_grad_hessian_ = 3.0;

  stan::test::expect_ad(tols, f(2), 4.0);
  stan::test::expect_ad(tols, f(4), 4.0);
  stan::test::expect_ad(tols, f(3), 5.0);

  stan::test::expect_ad(f(2), std::numeric_limits<double>::quiet_NaN());
}

TEST(mathMixScalFun, fallingFactorial_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::falling_factorial;
    return falling_factorial(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 0.5, 3.4;
  std::vector<int> std_in2{3, 1};
  stan::test::expect_ad_vectorized_binary(f, in1, std_in2);

  Eigen::MatrixXd mat_in1 = in1.replicate(1, 2);
  std::vector<std::vector<int>> std_std_in2{std_in2, std_in2};
  stan::test::expect_ad_vectorized_binary(f, mat_in1, std_std_in2);
}

TEST(mathMixScalFun, fallingFactorial_matvar) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::falling_factorial;
    return falling_factorial(x1, x2);
  };

  std::vector<int> std_in2{3, 1};
  Eigen::VectorXd in1(2);
  in1 << 0.5, 3.4;
  Eigen::MatrixXd mat(2, 2);
  mat << 0.5, 3.4, 0.5, 3.4;

  stan::test::expect_ad_matvar(f, in1, std_in2);
  stan::test::expect_ad_matvar(f, in1, std_in2[0]);
  stan::test::expect_ad_vectorized_matvar(f, mat, std_in2);
}
