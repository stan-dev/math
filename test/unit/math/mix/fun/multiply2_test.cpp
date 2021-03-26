#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_util.hpp>

TEST(mathMix, multiplicationPatterns2) {
  using stan::math::fvar;
  using stan::math::var;
  instantiate_multiply<fvar<fvar<double>>>();
  instantiate_multiply<fvar<var>>();
  instantiate_multiply<fvar<fvar<var>>>();
}
