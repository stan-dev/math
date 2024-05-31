#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cos) {
  auto f = [](const auto& x) {
    using stan::math::cos;
    return cos(x);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized<
      stan::test::ScalarSupport::RealAndComplex>(f, -2.6, -2, -0.2, -0.5, 0,
                                                 1.5, 3, 5, 5.3);
  stan::test::expect_complex_common(f);
}

TEST(mathMixMatFun, cos_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::cos;
    return cos(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-2.6, -2, -0.2, -0.5, 0, 1.5, 3, 5, 5.3};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
