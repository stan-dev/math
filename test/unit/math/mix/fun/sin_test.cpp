#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, sin) {
  auto f = [](const auto& x1) {
    using stan::math::sin;
    return sin(x1);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 3, 5, 5.3);
  stan::test::expect_complex_common(f);
}

TEST(mathMixMatFun, sin_varmat) {
  using stan::test::expect_ad_matvar;
  auto f = [](const auto& x1) {
    using stan::math::sin;
    return sin(x1);
  };
  Eigen::MatrixXd A(2, 3);
  A << -2.6, -2, -0.2, 3, 5, 5.3;
  expect_ad_matvar(f, A);
  std::vector<Eigen::MatrixXd> A_vec;
  A_vec.push_back(A);
  A_vec.push_back(A);
  A_vec.push_back(A);
  stan::test::expect_ad_matvar(f, A_vec);
}
