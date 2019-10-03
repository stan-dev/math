#include <stan/math/prim/arr.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>
#include <string>

/**
 * Require container
 */
TEST(requires, container_type_test) {
  EXPECT_FALSE((stan::container_value_type_check_base<
                stan::is_vector, std::is_floating_point, std::string>::value));
  EXPECT_TRUE(
      (stan::container_value_type_check_base<stan::is_vector,
                                             std::is_floating_point,
                                             std::vector<double>>::value));
}

template <template <class> class ContainerCheck,
          template <class> class TypeCheck, class Check, typename = void>
struct require_container_tester : std::false_type {};

template <template <class> class ContainerCheck, template <class> class TypeCheck, class Check>
struct require_container_tester<ContainerCheck, TypeCheck, Check, stan::require_container_vt<ContainerCheck, TypeCheck, Check>>
    : std::true_type {};

TEST(requires, generic_container_type_test) {
  EXPECT_FALSE(
      (require_container_tester<stan::is_vector, std::is_floating_point,
                                double>::value));
  EXPECT_TRUE((require_container_tester<stan::is_vector, std::is_floating_point,
                                        std::vector<double>>::value));
}


TEST(require_vt, std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt, std::vector>::unary<std::is_floating_point>();
}
TEST(require_vt, not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_vt, std::vector>::not_unary<std::is_floating_point>();
}
TEST(require_vt, all_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt, std::vector>::all<std::is_floating_point>();
}
TEST(require_vt, all_not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt, std::vector>::all_not<std::is_floating_point>();
}
TEST(require_vt, any_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt, std::vector>::any<std::is_floating_point>();
}
TEST(require_vt, any_not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt, std::vector>::any_not<std::is_floating_point>();
}

TEST(require_vt, vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_vt, std::vector>::unary<std::is_floating_point>();
}
TEST(require_vt, not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_vt, std::vector>::not_unary<std::is_floating_point>();
}
TEST(require_vt, all_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_vt, std::vector>::all<std::is_floating_point>();
}
TEST(require_vt, all_not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_vector_vt, std::vector>::all_not<std::is_floating_point>();
}
TEST(require_vt, any_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_vt, std::vector>::any<std::is_floating_point>();
}
TEST(require_vt, any_not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_vector_vt, std::vector>::any_not<std::is_floating_point>();
}

template <typename... Types>
using std_vector_vector = std::vector<std::vector<Types...>>;

TEST(require_container_vt, std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt, std_vector_vector>::unary<stan::is_std_vector>();
}
TEST(require_container_vt, not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_vt, std_vector_vector>::not_unary<stan::is_std_vector>();
}
TEST(require_container_vt, all_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt, std_vector_vector>::all<stan::is_std_vector>();
}
TEST(require_container_vt, all_not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt, std_vector_vector>::all_not<stan::is_std_vector>();
}
TEST(require_container_vt, any_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt, std_vector_vector>::any<stan::is_std_vector>();
}
TEST(require_container_vt, any_not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt, std_vector_vector>::any_not<stan::is_std_vector>();
}

TEST(require_container_vt, vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_vt, std_vector_vector>::unary<stan::is_vector>();
}
TEST(require_container_vt, not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_vt, std_vector_vector>::not_unary<stan::is_vector>();
}
TEST(require_container_vt, all_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_vt, std_vector_vector>::all<stan::is_vector>();
}
TEST(require_container_vt, all_not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_vector_vt, std_vector_vector>::all_not<stan::is_vector>();
}
TEST(require_container_vt, any_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_vt, std_vector_vector>::any<stan::is_vector>();
}
TEST(require_container_vt, any_not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_vector_vt, std_vector_vector>::any_not<stan::is_vector>();
}

TEST(require_st, std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_st, std_vector_vector>::unary<std::is_floating_point>();
}
TEST(require_st, not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_st, std_vector_vector>::not_unary<std::is_floating_point>();
}
TEST(require_st, all_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_st, std_vector_vector>::all<std::is_floating_point>();
}
TEST(require_st, all_not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_st, std_vector_vector>::all_not<std::is_floating_point>();
}
TEST(require_st, any_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_st, std_vector_vector>::any<std::is_floating_point>();
}
TEST(require_st, any_not_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_st, std_vector_vector>::any_not<std::is_floating_point>();
}

TEST(require_st, vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_st, std_vector_vector>::unary<std::is_floating_point>();
}
TEST(require_st, not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_st, std_vector_vector>::not_unary<std::is_floating_point>();
}
TEST(require_st, all_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_st, std_vector_vector>::all<std::is_floating_point>();
}
TEST(require_st, all_not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_vector_st, std_vector_vector>::all_not<std::is_floating_point>();
}
TEST(require_st, any_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_st, std_vector_vector>::any<std::is_floating_point>();
}
TEST(require_st, any_not_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_vector_st, std_vector_vector>::any_not<std::is_floating_point>();
}
