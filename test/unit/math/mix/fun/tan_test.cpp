#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, tan) {
  auto f = [](const auto& x) {
    using stan::math::tan;
    return tan(x);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2, -0.5, 0.5, 1.5, 3, 4.4);
  stan::test::expect_complex_common(f);
}

TEST(mathMixMatFun, sinh_varmat) {
  using stan::test::expect_ad_matvar;
  auto f = [](const auto& x1) {
    using stan::math::sinh;
    return sinh(x1);
  };
  Eigen::MatrixXd A(2, 3);
  A << -2, -0.5, 0.5, 1.5, 3, 4.4;
  expect_ad_matvar(f, A);
  std::vector<Eigen::MatrixXd> A_vec;
  A_vec.push_back(A);
  A_vec.push_back(A);
  A_vec.push_back(A);
  stan::test::expect_ad_matvar(f, A_vec);
}
