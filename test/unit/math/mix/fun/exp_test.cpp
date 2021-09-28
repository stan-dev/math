#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, exp) {
  auto f = [](const auto& x) {
    using stan::math::exp;
    return exp(x);
  };
  stan::test::expect_common_unary_vectorized<stan::test::PromoteToComplex::Yes>(
      f);
  stan::test::expect_unary_vectorized<stan::test::PromoteToComplex::Yes>(
      f, -15.2, -10, -0.5, 0.5, 1, 1.0, 1.3, 5, 10);
  stan::test::expect_complex_common(f);

  std::vector<double> com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> args{-2.6, -0.5, 0.5, 1.5};

  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(com_args));
  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(args));
}
