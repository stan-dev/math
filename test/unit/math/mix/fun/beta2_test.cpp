#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, beta_varmat_vectorized) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::beta;
    return beta(x1, x2);
  };

  Eigen::MatrixXd in1(2, 2);
  in1 << 0.5, 3.4, 5.2, 0.5;
  Eigen::MatrixXd in2(2, 2);
  in2 << 3.3, 0.9, 6.7, 3.3;
  stan::test::expect_ad_vectorized_matvar(f, in1, in2);
}
