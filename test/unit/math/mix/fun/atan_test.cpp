#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mathMixMatFun, atan) {
  auto f = [](const auto& x) {
    using stan::math::atan;
    return atan(x);
  };
  stan::test::expect_common_nonzero_unary_vectorized<
      stan::test::PromoteToComplex::No>(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 0.5, 1, 1.3, 1.5, 3);
  // avoid 0 imaginary component where autodiff doesn't work
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}

TEST(mathMixMatFun, atan_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::atan;
    return atan(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-2.6, -2, -0.2, 0.5, 1, 1.3, 1.5, 3};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
