#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <type_traits>

namespace stan {
namespace test {
namespace internal {
// taken from cppreference for is_detected
template <class T>
using copy_assign_t = decltype(std::declval<T&>() = std::declval<const T&>());

struct Meow {};
struct Purr {
  void operator=(const Purr&) = delete;
};

}  // namespace internal
}  // namespace test
}  // namespace stan

TEST(MetaTraits, is_detected_checks) {
  using stan::is_detected;
  using stan::test::internal::copy_assign_t;
  using stan::test::internal::Meow;
  using stan::test::internal::Purr;
  EXPECT_TRUE((is_detected<double, copy_assign_t>::value));
  EXPECT_TRUE((is_detected<int, copy_assign_t>::value));
  EXPECT_TRUE((is_detected<size_t, copy_assign_t>::value));
  EXPECT_TRUE((is_detected<Meow, copy_assign_t>::value));
  EXPECT_FALSE((is_detected<Purr, copy_assign_t>::value));
}
