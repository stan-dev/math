#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace stan {
namespace test {
template <typename... Types>
using check_fvar_or_arithmetic
    = math::conjunction<stan::is_fvar_or_arithmetic<Types>...>;
}
}  // namespace stan

TEST(MathMetaRevScal, is_fvar_or_arithmetic_simple) {
  using stan::is_fvar_or_arithmetic;
  EXPECT_TRUE(stan::is_fvar_or_arithmetic<stan::math::fvar<double>>::value);
  EXPECT_TRUE(stan::is_fvar_or_arithmetic<stan::math::fvar<double>&>::value);
  EXPECT_TRUE(stan::is_fvar_or_arithmetic<double>::value);
  EXPECT_TRUE(stan::is_fvar_or_arithmetic<double&>::value);
  bool temp = stan::test::check_fvar_or_arithmetic<
      stan::math::fvar<double>, std::vector<stan::math::fvar<double>>,
      stan::math::fvar<double>, stan::math::fvar<double>,
      stan::math::fvar<double>, stan::math::fvar<double>>::value;
  EXPECT_FALSE(temp);
  temp = stan::test::check_fvar_or_arithmetic<
      std::vector<stan::math::fvar<double>>, stan::math::fvar<double>,
      std::vector<stan::math::fvar<double> const*>>::value;
  EXPECT_FALSE(temp);
  temp = stan::test::check_fvar_or_arithmetic<
      std::vector<stan::math::fvar<double>>, stan::math::fvar<double>,
      std::vector<stan::math::fvar<double>>>::value;
  EXPECT_FALSE(temp);
}
