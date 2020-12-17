#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, tanh) {
  auto f = [](const auto& x1) {
    using stan::math::tanh;
    return tanh(x1);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1.2, -0.5, 0.5, 1.5);
  stan::test::expect_complex_common(f);
}

TEST(mathMixMatFun, tanh_varmat) {
  using stan::test::expect_ad_matvar;
  auto f = [](const auto& x1) {
    using stan::math::tanh;
    return tanh(x1);
  };
  Eigen::MatrixXd A(2, 3);
  A << -2.6, -2, -1.2, -0.5, 0.5, 1.5;
  expect_ad_matvar(f, A);
  std::vector<Eigen::MatrixXd> A_vec;
  A_vec.push_back(A);
  A_vec.push_back(A);
  A_vec.push_back(A);
  stan::test::expect_ad_matvar(f, A_vec);
}
