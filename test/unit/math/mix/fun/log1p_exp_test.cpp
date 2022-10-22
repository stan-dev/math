#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1pExp) {
  auto f = [](const auto& x1) { return stan::math::log1p_exp(x1); };
  stan::test::expect_common_nonzero_unary_vectorized<
      stan::test::ScalarSupport::Real>(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1, -0.5, -0.2, 0.5, 1.0,
                                      1.3, 2, 3);

  std::vector<double> com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> args{0.1, -2.5, 5.5};

  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(com_args));
  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(args));
}
