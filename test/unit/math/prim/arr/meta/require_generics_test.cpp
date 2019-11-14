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

template <template <class> class ContainerCheck,
          template <class> class TypeCheck, class Check>
struct require_container_tester<
    ContainerCheck, TypeCheck, Check,
    stan::require_container_vt<ContainerCheck, TypeCheck, Check>>
    : std::true_type {};

TEST(requires, generic_container_type_test) {
  EXPECT_FALSE(
      (require_container_tester<stan::is_vector, std::is_floating_point,
                                double>::value));
  EXPECT_TRUE((require_container_tester<stan::is_vector, std::is_floating_point,
                                        std::vector<double>>::value));
}

TEST(requires, std_vector_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_std_vector_t,
                       std::vector<double>>::unary();
  require_scal_checker<stan::require_not_std_vector_t, std::vector<double>,
                       std::vector<double>>::not_unary();
  require_scal_checker<stan::require_all_std_vector_t, std::vector<double>,
                       std::vector<double>>::all();
  require_scal_checker<stan::require_all_not_std_vector_t, std::vector<double>,
                       std::vector<double>>::all_not();
  require_scal_checker<stan::require_any_std_vector_t, std::vector<double>,
                       std::vector<double>>::any();
  require_scal_checker<stan::require_any_not_std_vector_t, std::vector<double>,
                       std::vector<double>>::any_not();
}

TEST(requires, vector_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_vector_t, std::vector<double>>::unary();
  require_scal_checker<stan::require_not_vector_t, std::vector<double>,
                       std::vector<double>>::not_unary();
  require_scal_checker<stan::require_all_vector_t, std::vector<double>,
                       std::vector<double>>::all();
  require_scal_checker<stan::require_all_not_vector_t, std::vector<double>,
                       std::vector<double>>::all_not();
  require_scal_checker<stan::require_any_vector_t, std::vector<double>,
                       std::vector<double>>::any();
  require_scal_checker<stan::require_any_not_vector_t, std::vector<double>,
                       std::vector<double>>::any_not();
}

TEST(requires, vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_vector_like_t,
                       std::vector<double>>::unary();
  require_scal_checker<stan::require_not_vector_like_t, std::vector<double>,
                       std::vector<double>>::not_unary();
  require_scal_checker<stan::require_all_vector_like_t, std::vector<double>,
                       std::vector<double>>::all();
  require_scal_checker<stan::require_all_not_vector_like_t, std::vector<double>,
                       std::vector<double>>::all_not();
  require_scal_checker<stan::require_any_vector_like_t, std::vector<double>,
                       std::vector<double>>::any();
  require_scal_checker<stan::require_any_not_vector_like_t, std::vector<double>,
                       std::vector<double>>::any_not();
}

TEST(requires, std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt,
                            std::vector>::unary<std::is_floating_point>();
}
TEST(requires, not_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_vt,
                            std::vector>::not_unary<std::is_floating_point>();
}
TEST(requires, all_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt,
                            std::vector>::all<std::is_floating_point>();
}
TEST(requires, all_not_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt,
                            std::vector>::all_not<std::is_floating_point>();
}
TEST(requires, any_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt,
                            std::vector>::any<std::is_floating_point>();
}
TEST(requires, any_not_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt,
                            std::vector>::any_not<std::is_floating_point>();
}

TEST(requires, vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_vt,
                            std::vector>::unary<std::is_floating_point>();
}
TEST(requires, not_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_vt,
                            std::vector>::not_unary<std::is_floating_point>();
}
TEST(requires, all_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_vt,
                            std::vector>::all<std::is_floating_point>();
}
TEST(requires, all_not_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_vector_vt,
                            std::vector>::all_not<std::is_floating_point>();
}
TEST(requires, any_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_vt,
                            std::vector>::any<std::is_floating_point>();
}
TEST(requires, any_not_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_vector_vt,
                            std::vector>::any_not<std::is_floating_point>();
}

template <typename... Types>
using std_vector_vector = std::vector<std::vector<Types...>>;

TEST(requires, std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt,
                            std_vector_vector>::unary<stan::is_std_vector>();
}
TEST(requires, not_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_std_vector_vt,
      std_vector_vector>::not_unary<stan::is_std_vector>();
}
TEST(requires, all_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt,
                            std_vector_vector>::all<stan::is_std_vector>();
}
TEST(requires, all_not_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt,
                            std_vector_vector>::all_not<stan::is_std_vector>();
}
TEST(requires, any_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt,
                            std_vector_vector>::any<stan::is_std_vector>();
}
TEST(requires, any_not_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt,
                            std_vector_vector>::any_not<stan::is_std_vector>();
}

TEST(requires, vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_vt,
                            std_vector_vector>::unary<stan::is_vector>();
}
TEST(requires, not_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_vt,
                            std_vector_vector>::not_unary<stan::is_vector>();
}
TEST(requires, all_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_vt,
                            std_vector_vector>::all<stan::is_vector>();
}
TEST(requires, all_not_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_vector_vt,
                            std_vector_vector>::all_not<stan::is_vector>();
}
TEST(requires, any_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_vt,
                            std_vector_vector>::any<stan::is_vector>();
}
TEST(requires, any_not_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_vector_vt,
                            std_vector_vector>::any_not<stan::is_vector>();
}

TEST(requires, std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_st,
                            std_vector_vector>::unary<std::is_floating_point>();
}
TEST(requires, not_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_std_vector_st,
      std_vector_vector>::not_unary<std::is_floating_point>();
}
TEST(requires, all_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_st,
                            std_vector_vector>::all<std::is_floating_point>();
}
TEST(requires, all_not_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_std_vector_st,
      std_vector_vector>::all_not<std::is_floating_point>();
}
TEST(requires, any_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_st,
                            std_vector_vector>::any<std::is_floating_point>();
}
TEST(requires, any_not_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_std_vector_st,
      std_vector_vector>::any_not<std::is_floating_point>();
}

TEST(requires, vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_st,
                            std_vector_vector>::unary<std::is_floating_point>();
}
TEST(requires, not_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_st, std_vector_vector>::
      not_unary<std::is_floating_point>();
}
TEST(requires, all_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_st,
                            std_vector_vector>::all<std::is_floating_point>();
}
TEST(requires, all_not_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_vector_st,
      std_vector_vector>::all_not<std::is_floating_point>();
}
TEST(requires, any_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_st,
                            std_vector_vector>::any<std::is_floating_point>();
}
TEST(requires, any_not_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_vector_st,
      std_vector_vector>::any_not<std::is_floating_point>();
}
