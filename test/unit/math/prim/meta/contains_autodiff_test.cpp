#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

TEST(MathMetaContainsAutodiff, scalar_detection) {
  using stan::math::fvar;
  EXPECT_FALSE(stan::contains_autodiff_v<double>);
  EXPECT_FALSE(stan::contains_autodiff_v<int>);
  EXPECT_TRUE(stan::contains_autodiff_v<fvar<double>>);
}

TEST(MathMetaContainsAutodiff, container_detection) {
  using stan::math::fvar;
  EXPECT_FALSE((stan::contains_autodiff_v<std::vector<double>>));
  EXPECT_TRUE((stan::contains_autodiff_v<std::vector<fvar<double>>>));
  EXPECT_FALSE(
      (stan::contains_autodiff_v<std::tuple<int, std::vector<double>>>));
  EXPECT_TRUE((stan::contains_autodiff_v<
               std::tuple<int, std::vector<fvar<double>>>>));
}

TEST(MathMetaContainsAutodiff, nested_tuple_detection) {
  using stan::math::fvar;
  using nested_ad_t = std::tuple<double, std::tuple<int, fvar<double>>>;
  using nested_non_ad_t = std::tuple<double, std::tuple<int, double>>;
  EXPECT_TRUE((stan::contains_autodiff_v<nested_ad_t>));
  EXPECT_FALSE((stan::contains_autodiff_v<nested_non_ad_t>));
}
