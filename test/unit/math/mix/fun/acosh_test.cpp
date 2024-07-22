#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, acosh) {
  auto f = [](const auto& x1) {
    using stan::math::acosh;
    return acosh(x1);
  };
  for (double x : stan::test::internal::common_args())
    stan::test::expect_unary_vectorized(f, x);
  stan::test::expect_unary_vectorized<
      stan::test::ScalarSupport::RealAndComplex>(f, 1.5, 3.2, 5, 10, 12.9);
  // avoid pole at complex zero that can't be autodiffed
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}

TEST(mathMixMatFun, acosh_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::acosh;
    return acosh(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{1.5, 3.2, 5, 10, 12.9};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
