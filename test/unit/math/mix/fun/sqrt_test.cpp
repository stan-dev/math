#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, sqrt) {
  auto f = [](const auto& x1) {
    using stan::math::sqrt;
    return sqrt(x1);
  };
  stan::test::expect_common_nonzero_unary_vectorized<
      stan::test::ScalarSupport::Real>(f);
  stan::test::expect_unary_vectorized(f, -6, -5.2, 1.3, 7, 10.7, 36, 1e6);

  // undefined with 0 in denominator
  stan::test::expect_ad(f, std::complex<double>(0.9, 0.8));
  for (double im : std::vector<double>{-1.3, 2.3}) {
    for (double re : std::vector<double>{-3.6, -0.0, 0.0, 0.5}) {
      stan::test::expect_ad(f, std::complex<double>(re, im));
    }
  }
}

TEST(mathMixMatFun, sqrt_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::sqrt;
    return sqrt(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-6, -5.2, 1.3, 7, 10.7, 36, 1e6};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
