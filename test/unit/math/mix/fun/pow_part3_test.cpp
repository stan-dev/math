#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

TEST(mathMixScalFun, pow_varmat) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::pow;
    using std::pow;
    return pow(x1, x2);
  };
  Eigen::MatrixXd mat1(2, 4);
  mat1 << -0.4, 0.5, 0.5, 0.5, 0.5, 1.0, 3.0, 4.0;
  Eigen::MatrixXd mat2(2, 4);
  mat2 << 0.5, 0.5, 1.0, 1.2, 5.0, 2.0, 4.0, -2.0;
  stan::test::expect_ad_matvar(f, mat1, mat2);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad_matvar(f, mat1, nan);
  stan::test::expect_ad_matvar(f, nan, mat2);

  Eigen::VectorXd in1(3);
  in1 << 0.5, 3.4, 5.2;
  Eigen::VectorXd in2(3);
  in2 << 3.3, 0.9, 2.1;
  stan::test::expect_ad_matvar(f, in1, in2);
  stan::test::expect_ad_matvar(f, in1, 2.0);
  stan::test::expect_ad_matvar(f, 2.0, in1);

  Eigen::MatrixXd mat_in1(2, 2);
  mat_in1 << 0.5, 3.4, 0.5, 3.4;
  std::vector<int> std_in2{3, 1};
  stan::test::expect_ad_vectorized_matvar(f, mat_in1, std_in2);
}
