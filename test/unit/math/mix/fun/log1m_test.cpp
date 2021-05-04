#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1m) {
  auto f = [](const auto& x1) { return stan::math::log1m(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -21.5, -21, -1, 0.9, 3);

  std::vector<double> com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> args{-2.6, -0.5, 0.5, 1.5};

  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(com_args));
  stan::test::expect_ad_vector_matvar(f, stan::math::to_vector(args));
}
