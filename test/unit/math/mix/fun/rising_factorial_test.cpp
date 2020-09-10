#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, risingFactorial) {
  auto f = [](const int x2) {
    return [=](const auto& x1) { return stan::math::rising_factorial(x1, x2); };
  };

  stan::test::expect_ad(f(1), 4.0);
  stan::test::expect_ad(f(3), 5.0);
  stan::test::expect_ad(f(4), 4.0);

  // 3rd derivatives close to zero and rel tolerance fails
  stan::test::ad_tolerances tols;
  tols.grad_hessian_grad_hessian_ = 3.0;

  stan::test::expect_ad(tols, f(1), 5.0);
  stan::test::expect_ad(tols, f(2), 4.0);
  stan::test::expect_ad(tols, f(2), 4.0);
  stan::test::expect_ad(tols, f(4), 4.0);

  stan::test::expect_ad(f(2), std::numeric_limits<double>::quiet_NaN());
}

TEST(mathMixScalFun, risingFactorial_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::rising_factorial;
    return rising_factorial(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 0.5, 3.4;
  std::vector<int> std_in2{3, 1};
  stan::test::expect_ad_vectorized_binary(f, in1, std_in2);

  Eigen::MatrixXd mat_in1 = in1.replicate(1, 2);
  std::vector<std::vector<int>> std_std_in2{std_in2, std_in2};
  stan::test::expect_ad_vectorized_binary(f, mat_in1, std_std_in2);
}
