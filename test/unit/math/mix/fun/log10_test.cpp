#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, log10) {
  auto f = [](const auto& x1) {
    using stan::math::log10;
    return log10(x1);
  };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -0.2, 1e-3, 1, 1.3, 3, 3.7, 10, 10.2,
                                      1e6);

  // non-zero real and imaginary components
  for (auto re : std::vector<double>{-2.7, 1, 2.3}) {
    for (auto im : std::vector<double>{-1.5, 1.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
  // zero tests which lead to finite vals
  stan::test::expect_ad(f, std::complex<double>{0, 2.1});
  stan::test::expect_ad(f, std::complex<double>{0, -2.1});
  stan::test::expect_ad(f, std::complex<double>{2.1, 0});
  stan::test::expect_ad(f, std::complex<double>{-0.0, 2.1});
  stan::test::expect_ad(f, std::complex<double>{-0.0, -2.1});
  stan::test::expect_ad(f, std::complex<double>{2.1, -0.0});
  // (negative real and zero imaginary illegal)
}

TEST(mathMixMatFun, log10_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::log10;
    return log10(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-0.2, 1e-3, 1, 1.3, 3, 3.7, 10, 10.2, 1e6};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
