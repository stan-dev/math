#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace stan {
namespace math {
namespace test {
template <typename T>
class dummy_base {};
class dummy_derived : public dummy_base<dummy_derived> {};
}  // namespace test
}  // namespace math
}  // namespace stan

TEST(MathMetaPrim, is_convertible_to_checks) {
  using stan::is_base_pointer_convertible;
  using stan::math::test::dummy_base;
  using stan::math::test::dummy_derived;
  EXPECT_TRUE((is_base_pointer_convertible<dummy_base, dummy_derived>::value));
}
