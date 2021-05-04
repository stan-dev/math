#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <type_traits>

namespace stan {
namespace test {
namespace internal {

// taken from cppreference for void_t
template <typename T, typename = void>
struct is_iterable : std::false_type {};
template <typename T>
struct is_iterable<T, stan::void_t<decltype(std::declval<T>().begin()),
                                   decltype(std::declval<T>().end())>>
    : std::true_type {};

class dummy_iterable {
 public:
  int begin() { return 1; }
  int end() { return 1; }
};
}  // namespace internal
}  // namespace test
}  // namespace stan

TEST(MetaTraits, void_t_checks) {
  using stan::void_t;
  using stan::test::internal::dummy_iterable;
  using stan::test::internal::is_iterable;
  EXPECT_FALSE(is_iterable<double>::value);
  EXPECT_FALSE(is_iterable<int>::value);
  EXPECT_FALSE(is_iterable<size_t>::value);
  EXPECT_TRUE(is_iterable<dummy_iterable>::value);
}
