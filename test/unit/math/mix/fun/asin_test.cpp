#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, asin) {
  auto f = [](const auto& x) {
    using stan::math::asin;
    return asin(x);
  };
  // can't autodiff asin through integers
  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 0.7, 1.3, 3.4, 5);
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}

TEST(mathMixMatFun, asin_varmat) {
  using stan::test::expect_ad_matvar;
  auto f = [](const auto& x1) {
    using stan::math::asin;
    return asin(x1);
  };
  Eigen::VectorXd A(7);
  A << -2.6, -2, -0.2, 0.7, 1.3, 3.4, 5;
  expect_ad_matvar(f, A);
}
