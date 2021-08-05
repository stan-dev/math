#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1p) {
  auto f = [](const auto& x1) { return stan::math::log1p(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2, -.99, -0.5, 0.5, 7.2, 1000);

  std::vector<double> com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> args{0.1, -0.5, 5.5};
  std::vector<double> args_invalid{-2.5, -1.1};

  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(com_args));
  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(args));
  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(args_invalid));
}
