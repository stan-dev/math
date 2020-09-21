#include <stan/math/rev/meta.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>
#include <string>

// Basic tests for the underlying requires
TEST(requires_prim_scal, t_test) {
  using stan::require_t;
  using stan::test::unary_require_tester;
  EXPECT_TRUE((unary_require_tester<require_t, std::true_type>::value));
  EXPECT_FALSE((unary_require_tester<require_t, std::false_type>::value));
}

TEST(requires_prim_scal, not_test) {
  using stan::require_not_t;
  using stan::test::unary_require_tester;
  EXPECT_TRUE((unary_require_tester<require_not_t, std::false_type>::value));
  EXPECT_FALSE((unary_require_tester<require_not_t, std::true_type>::value));
}

TEST(requires_prim_scal, all_not_test) {
  using stan::require_all_not_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_all_not_t, std::false_type,
                                        std::false_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_all_not_t, std::true_type,
                                         std::true_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_all_not_t, std::true_type,
                                         std::false_type>::value));
}

TEST(requires_prim_scal, any_not_test) {
  using stan::require_any_not_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_any_not_t, std::false_type,
                                        std::false_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_t, std::true_type,
                                         std::true_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_any_not_t, std::true_type,
                                        std::false_type>::value));
}

TEST(requires_prim_scal, all_test) {
  using stan::require_all_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_t, std::false_type,
                                         std::false_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_t, std::true_type,
                                        std::true_type>::value));
  EXPECT_FALSE((require_variadic_checker<require_all_t, std::true_type,
                                         std::false_type>::value));
}

TEST(requires_prim_scal, any_test) {
  using stan::require_any_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_any_t, std::false_type,
                                         std::false_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_any_t, std::true_type,
                                        std::true_type>::value));
  EXPECT_TRUE((require_variadic_checker<require_any_t, std::true_type,
                                        std::false_type>::value));
}

// Test same
TEST(requires_prim_scal, same_test) {
  using stan::require_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE(
      (require_variadic_checker<stan::require_same_t, double, int>::value));
  EXPECT_TRUE(
      (require_variadic_checker<stan::require_same_t, double, double>::value));
  EXPECT_TRUE(
      (require_variadic_checker<stan::require_same_t, int, int>::value));
}

TEST(requires_prim_scal, not_same_test) {
  using stan::require_not_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE(
      (require_variadic_checker<require_not_same_t, double, int>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_not_same_t, double, double>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_same_t, int, int>::value));
}

TEST(requires_prim_scal, all_same_test) {
  using stan::require_all_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_same_t, double,
                                         std::string, double>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_same_t, double, double,
                                        double>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_all_same_t, int, int, int>::value));
}

TEST(requires_prim_scal, any_not_same_test) {
  using stan::require_any_not_same_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_any_not_same_t, double, int,
                                        int>::value));
  EXPECT_TRUE((require_variadic_checker<require_any_not_same_t, double, int,
                                        double>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_same_t, double, double,
                                         double>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_any_not_same_t, int, int, int>::value));
}

// Test convertible
TEST(requires_prim_scal, convertible_test) {
  using stan::require_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<stan::require_convertible_t, double,
                                         char[1]>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_convertible_t, double,
                                        double>::value));
  EXPECT_TRUE(
      (require_variadic_checker<stan::require_convertible_t, int, int>::value));
}

TEST(requires_prim_scal, not_convertible_test) {
  using stan::require_not_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_not_convertible_t, double,
                                        char[1]>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_convertible_t, double,
                                         double>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_not_convertible_t, int, int>::value));
}

TEST(requires_prim_scal, all_convertible_test) {
  using stan::require_all_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_convertible_t, double,
                                         std::string, std::string>::value));
  EXPECT_FALSE((require_variadic_checker<require_all_convertible_t, double,
                                         std::string, double>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_convertible_t, double, int,
                                        int>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_convertible_t, int, int,
                                        int>::value));
}

TEST(requires_prim_scal, any_not_convertible_test) {
  using stan::require_any_not_convertible_t;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_any_not_convertible_t, double,
                                        std::string, std::string>::value));
  EXPECT_TRUE((require_variadic_checker<require_any_not_convertible_t, double,
                                        int, std::string>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_convertible_t, double,
                                         double, double>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_convertible_t, int,
                                         double, int>::value));
}

// Double or Int
TEST(requires_prim_scal, double_or_int_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_double_or_int_t, double, int>::unary();
  require_scal_checker<stan::require_not_double_or_int_t, double,
                       int>::not_unary();
  require_scal_checker<stan::require_all_double_or_int_t, double, int>::all();
  require_scal_checker<stan::require_all_not_double_or_int_t, double,
                       int>::all_not();
  require_scal_checker<stan::require_any_double_or_int_t, double, int>::any();
  require_scal_checker<stan::require_any_not_double_or_int_t, double,
                       int>::any_not();
}

// Double or Int
TEST(requires_prim_scal, require_string_convertible_test) {
  using stan::test::require_scal_checker;
  using std::string;
  require_scal_checker<stan::require_string_convertible_t, string,
                       const char*>::unary();
  require_scal_checker<stan::require_not_string_convertible_t, string,
                       const char*>::not_unary();
  require_scal_checker<stan::require_all_string_convertible_t, string,
                       const char*>::all();
  require_scal_checker<stan::require_all_not_string_convertible_t, string,
                       const char*>::all_not();
  require_scal_checker<stan::require_any_string_convertible_t, string,
                       const char*>::any();
  require_scal_checker<stan::require_any_not_string_convertible_t, string,
                       const char*>::any_not();
}

TEST(requires_prim_scal, arithmetic_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_arithmetic_t, double, float, int>::unary();
  require_scal_checker<stan::require_not_arithmetic_t, double, float,
                       int>::not_unary();
  require_scal_checker<stan::require_all_arithmetic_t, double, float,
                       int>::all();
  require_scal_checker<stan::require_all_not_arithmetic_t, double, float,
                       int>::all_not();
  require_scal_checker<stan::require_any_arithmetic_t, double, float,
                       int>::any();
  require_scal_checker<stan::require_any_not_arithmetic_t, double, float,
                       int>::any_not();
}

TEST(requires_prim_scal, var_or_arithmetic_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_var_or_arithmetic_t, double, int>::unary();
  require_scal_checker<stan::require_not_var_or_arithmetic_t, double,
                       int>::not_unary();
  require_scal_checker<stan::require_all_var_or_arithmetic_t, double,
                       int>::all();
  require_scal_checker<stan::require_all_not_var_or_arithmetic_t, double,
                       int>::all_not();
  require_scal_checker<stan::require_any_var_or_arithmetic_t, double,
                       int>::any();
  require_scal_checker<stan::require_any_not_var_or_arithmetic_t, double,
                       int>::any_not();
}

TEST(requires_prim_arr, std_vector_t_test) {
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

TEST(requires_prim_arr, vector_t_test) {
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

TEST(requires_prim_arr, vector_like_t_test) {
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

TEST(requires_prim_arr, std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt,
                            std::vector>::unary<std::is_floating_point>();
}
TEST(requires_prim_arr, not_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_vt,
                            std::vector>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_arr, all_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt,
                            std::vector>::all<std::is_floating_point>();
}
TEST(requires_prim_arr, all_not_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt,
                            std::vector>::all_not<std::is_floating_point>();
}
TEST(requires_prim_arr, any_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt,
                            std::vector>::any<std::is_floating_point>();
}
TEST(requires_prim_arr, any_not_std_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt,
                            std::vector>::any_not<std::is_floating_point>();
}

TEST(requires_prim_arr, vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_vt,
                            std::vector>::unary<std::is_floating_point>();
}
TEST(requires_prim_arr, not_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_vt,
                            std::vector>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_arr, all_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_vt,
                            std::vector>::all<std::is_floating_point>();
}
TEST(requires_prim_arr, all_not_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_vector_vt,
                            std::vector>::all_not<std::is_floating_point>();
}
TEST(requires_prim_arr, any_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_vt,
                            std::vector>::any<std::is_floating_point>();
}
TEST(requires_prim_arr, any_not_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_vector_vt,
                            std::vector>::any_not<std::is_floating_point>();
}

template <typename... Types>
using std_vector_vector = std::vector<std::vector<Types...>>;

TEST(requires_prim_arr, std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt,
                            std_vector_vector>::unary<stan::is_std_vector>();
}
TEST(requires_prim_arr, not_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_std_vector_vt,
      std_vector_vector>::not_unary<stan::is_std_vector>();
}
TEST(requires_prim_arr, all_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt,
                            std_vector_vector>::all<stan::is_std_vector>();
}
TEST(requires_prim_arr, all_not_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt,
                            std_vector_vector>::all_not<stan::is_std_vector>();
}
TEST(requires_prim_arr, any_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt,
                            std_vector_vector>::any<stan::is_std_vector>();
}
TEST(requires_prim_arr, any_not_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt,
                            std_vector_vector>::any_not<stan::is_std_vector>();
}

TEST(requires_prim_arr, vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_vt,
                            std_vector_vector>::unary<stan::is_vector>();
}
TEST(requires_prim_arr, not_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_vt,
                            std_vector_vector>::not_unary<stan::is_vector>();
}
TEST(requires_prim_arr, all_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_vt,
                            std_vector_vector>::all<stan::is_vector>();
}
TEST(requires_prim_arr, all_not_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_vector_vt,
                            std_vector_vector>::all_not<stan::is_vector>();
}
TEST(requires_prim_arr, any_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_vt,
                            std_vector_vector>::any<stan::is_vector>();
}
TEST(requires_prim_arr, any_not_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_vector_vt,
                            std_vector_vector>::any_not<stan::is_vector>();
}

TEST(requires_prim_arr, std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_st,
                            std_vector_vector>::unary<std::is_floating_point>();
}
TEST(requires_prim_arr, not_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_std_vector_st,
      std_vector_vector>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_arr, all_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_st,
                            std_vector_vector>::all<std::is_floating_point>();
}
TEST(requires_prim_arr, all_not_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_std_vector_st,
      std_vector_vector>::all_not<std::is_floating_point>();
}
TEST(requires_prim_arr, any_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_st,
                            std_vector_vector>::any<std::is_floating_point>();
}
TEST(requires_prim_arr, any_not_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_std_vector_st,
      std_vector_vector>::any_not<std::is_floating_point>();
}

TEST(requires_prim_arr, vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_vector_st,
                            std_vector_vector>::unary<std::is_floating_point>();
}
TEST(requires_prim_arr, not_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_vector_st, std_vector_vector>::
      not_unary<std::is_floating_point>();
}
TEST(requires_prim_arr, all_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_vector_st,
                            std_vector_vector>::all<std::is_floating_point>();
}
TEST(requires_prim_arr, all_not_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_vector_st,
      std_vector_vector>::all_not<std::is_floating_point>();
}
TEST(requires_prim_arr, any_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_vector_st,
                            std_vector_vector>::any<std::is_floating_point>();
}
TEST(requires_prim_arr, any_not_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_vector_st,
      std_vector_vector>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_x = Eigen::Matrix<T, -1, -1>;

// Test same
TEST(requires_prim_mat, same_vtest) {
  using stan::require_vt_same;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<stan::require_st_same, eigen_x<double>,
                                         eigen_x<int>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_st_same, eigen_x<double>,
                                        eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_st_same, eigen_x<int>,
                                        eigen_x<int>>::value));
}

TEST(requires_prim_mat, not_same_vtest) {
  using stan::require_not_st_same;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_not_st_same, eigen_x<double>,
                                        eigen_x<int>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_st_same, eigen_x<double>,
                                         eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_st_same, eigen_x<int>,
                                         eigen_x<int>>::value));
}

TEST(requires_prim_mat, all_same_vtest) {
  using stan::require_all_st_same;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_st_same, eigen_x<double>,
                                         std::string, eigen_x<double>>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_all_st_same, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_st_same, eigen_x<int>,
                                        eigen_x<int>, eigen_x<int>>::value));
}

TEST(requires_prim_mat, any_not_same_vtest) {
  using stan::require_any_not_st_same;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE(
      (require_variadic_checker<require_any_not_st_same, eigen_x<double>,
                                eigen_x<int>, eigen_x<int>>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_any_not_st_same, eigen_x<double>,
                                eigen_x<int>, eigen_x<double>>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_any_not_st_same, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_st_same, eigen_x<int>,
                                         eigen_x<int>, eigen_x<int>>::value));
}

// Test same
TEST(requires_prim_mat, same_stest) {
  using stan::require_st_same;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<stan::require_st_same, eigen_x<double>,
                                         eigen_x<int>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_st_same, eigen_x<double>,
                                        eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_st_same, eigen_x<int>,
                                        eigen_x<int>>::value));
}

TEST(requires_prim_mat, not_same_stest) {
  using stan::require_not_st_same;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_not_st_same, eigen_x<double>,
                                        eigen_x<int>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_st_same, eigen_x<double>,
                                         eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_st_same, eigen_x<int>,
                                         eigen_x<int>>::value));
}

TEST(requires_prim_mat, all_same_stest) {
  using stan::require_all_st_same;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_st_same, eigen_x<double>,
                                         std::string, eigen_x<double>>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_all_st_same, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_st_same, eigen_x<int>,
                                        eigen_x<int>, eigen_x<int>>::value));
}

TEST(requires_prim_mat, any_not_same_stest) {
  using stan::require_any_not_st_same;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE(
      (require_variadic_checker<require_any_not_st_same, eigen_x<double>,
                                eigen_x<int>, eigen_x<int>>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_any_not_st_same, eigen_x<double>,
                                eigen_x<int>, eigen_x<double>>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_any_not_st_same, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_st_same, eigen_x<int>,
                                         eigen_x<int>, eigen_x<int>>::value));
}

TEST(requires_prim_mat, eigen_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_eigen_t, eigen_x<double>>::unary();
}
TEST(requires_prim_mat, eigen_not_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_eigen_t, eigen_x<double>,
                       eigen_x<double>>::not_unary();
}
TEST(requires_prim_mat, eigen_all_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_eigen_t, eigen_x<double>,
                       eigen_x<double>>::all();
}
TEST(requires_prim_mat, eigen_all_not_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_eigen_t, eigen_x<double>,
                       eigen_x<double>>::all_not();
}
TEST(requires_prim_mat, eigen_any_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_eigen_t, eigen_x<double>,
                       eigen_x<double>>::any();
}
TEST(requires_prim_mat, eigen_any_not_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_eigen_t, eigen_x<double>,
                       eigen_x<double>>::any_not();
}

TEST(requires_prim_mat, eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vt,
                            eigen_x>::unary<std::is_floating_point>();
}
TEST(requires_prim_mat, not_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vt,
                            eigen_x>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_mat, all_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vt,
                            eigen_x>::all<std::is_floating_point>();
}
TEST(requires_prim_mat, all_not_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vt,
                            eigen_x>::all_not<std::is_floating_point>();
}
TEST(requires_prim_mat, any_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vt,
                            eigen_x>::any<std::is_floating_point>();
}
TEST(requires_prim_mat, any_not_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vt,
                            eigen_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_vector_x = Eigen::Matrix<T, 1, -1>;

TEST(requires_prim_mat, eigen_vector_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_eigen_vector_t,
                       eigen_vector_x<double>>::unary();
  require_scal_checker<stan::require_not_eigen_vector_t, eigen_vector_x<double>,
                       eigen_vector_x<double>>::not_unary();
  require_scal_checker<stan::require_all_eigen_vector_t, eigen_vector_x<double>,
                       eigen_vector_x<double>>::all();
  require_scal_checker<stan::require_all_not_eigen_vector_t,
                       eigen_vector_x<double>,
                       eigen_vector_x<double>>::all_not();
  require_scal_checker<stan::require_any_eigen_vector_t, eigen_vector_x<double>,
                       eigen_vector_x<double>>::any();
  require_scal_checker<stan::require_any_not_eigen_vector_t,
                       eigen_vector_x<double>,
                       eigen_vector_x<double>>::any_not();
}

TEST(requires_prim_mat, eigen_std_vector_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_vector_t, std::vector<double>,
                       eigen_vector_x<double>>::unary();
  require_scal_checker<stan::require_not_vector_t, std::vector<double>,
                       eigen_vector_x<double>>::not_unary();
  require_scal_checker<stan::require_all_vector_t, std::vector<double>,
                       eigen_vector_x<double>>::all();
  require_scal_checker<stan::require_all_not_vector_t, std::vector<double>,
                       eigen_vector_x<double>>::all_not();
  require_scal_checker<stan::require_any_vector_t, std::vector<double>,
                       eigen_vector_x<double>>::any();
  require_scal_checker<stan::require_any_not_vector_t, std::vector<double>,
                       eigen_vector_x<double>>::any_not();
}

TEST(requires_prim_mat, eigen_vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_vector_like_t, std::vector<double>,
                       eigen_vector_x<double>>::unary();
}

TEST(requires, not_eigen_vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_vector_like_t, std::vector<double>,
                       eigen_vector_x<double>>::not_unary();
}
TEST(requires, all_eigen_vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_vector_like_t, std::vector<double>,
                       eigen_vector_x<double>>::all();
}
TEST(requires, all_not_eigen_vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_vector_like_t, std::vector<double>,
                       eigen_vector_x<double>>::all_not();
}
TEST(requires, any_eigen_vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_vector_like_t, std::vector<double>,
                       eigen_vector_x<double>>::any();
}
TEST(requires, any_not_eigen_vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_vector_like_t, std::vector<double>,
                       std::vector<double>>::any_not();
}

TEST(requires_prim_mat, eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_vector_x>::unary<std::is_floating_point>();
}
TEST(requires_prim_mat, not_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vector_vt, eigen_vector_x>::
      not_unary<std::is_floating_point>();
}
TEST(requires_prim_mat, all_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_vector_x>::all<std::is_floating_point>();
}
TEST(requires_prim_mat, all_not_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_vector_x>::all_not<std::is_floating_point>();
}
TEST(requires_prim_mat, any_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_vector_x>::any<std::is_floating_point>();
}
TEST(requires_prim_mat, any_not_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_eigen_x = Eigen::Matrix<Eigen::Matrix<T, -1, -1>, -1, -1>;

TEST(requires_prim_mat, eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vt,
                            eigen_eigen_x>::unary<stan::is_eigen>();
}
TEST(requires_prim_mat, not_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vt,
                            eigen_eigen_x>::not_unary<stan::is_eigen>();
}
TEST(requires_prim_mat, all_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vt,
                            eigen_eigen_x>::all<stan::is_eigen>();
}
TEST(requires_prim_mat, all_not_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vt,
                            eigen_eigen_x>::all_not<stan::is_eigen>();
}
TEST(requires_prim_mat, any_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vt,
                            eigen_eigen_x>::any<stan::is_eigen>();
}
TEST(requires_prim_mat, any_not_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vt,
                            eigen_eigen_x>::any_not<stan::is_eigen>();
}

TEST(requires_prim_mat, eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_st,
                            eigen_eigen_x>::unary<std::is_floating_point>();
}
TEST(requires_prim_mat, not_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_st,
                            eigen_eigen_x>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_mat, all_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_st,
                            eigen_eigen_x>::all<std::is_floating_point>();
}
TEST(requires_prim_mat, all_not_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_st,
                            eigen_eigen_x>::all_not<std::is_floating_point>();
}
TEST(requires_prim_mat, any_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_st,
                            eigen_eigen_x>::any<std::is_floating_point>();
}
TEST(requires_prim_mat, any_not_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_st,
                            eigen_eigen_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_eigen_vector_x = Eigen::Matrix<Eigen::Matrix<T, 1, -1>, 1, -1>;

TEST(requires_prim_mat, eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_eigen_vector_x>::unary<stan::is_vector>();
}
TEST(requires_prim_mat, not_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::not_unary<stan::is_vector>();
}
TEST(requires_prim_mat, all_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_eigen_vector_x>::all<stan::is_vector>();
}
TEST(requires_prim_mat, all_not_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::all_not<stan::is_vector>();
}
TEST(requires_prim_mat, any_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_eigen_vector_x>::any<stan::is_vector>();
}
TEST(requires_prim_mat, any_not_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::any_not<stan::is_vector>();
}

TEST(requires_prim_mat, eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_eigen_vector_st,
      eigen_eigen_vector_x>::unary<std::is_floating_point>();
}
TEST(requires_prim_mat, not_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_st,
      eigen_eigen_vector_x>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_mat, all_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_eigen_vector_st,
      eigen_eigen_vector_x>::all<std::is_floating_point>();
}
TEST(requires_prim_mat, all_not_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_eigen_vector_st,
      eigen_eigen_vector_x>::all_not<std::is_floating_point>();
}
TEST(requires_prim_mat, any_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_eigen_vector_st,
      eigen_eigen_vector_x>::any<std::is_floating_point>();
}
TEST(requires_prim_mat, any_not_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_eigen_vector_st,
      eigen_eigen_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_std_vector_x = Eigen::Matrix<std::vector<T>, 1, -1>;

TEST(requires_prim_mat, eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_std_vector_x>::unary<stan::is_std_vector>();
}
TEST(requires_prim_mat, not_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_vt,
      eigen_std_vector_x>::not_unary<stan::is_std_vector>();
}
TEST(requires_prim_mat, all_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_std_vector_x>::all<stan::is_std_vector>();
}
TEST(requires_prim_mat, all_not_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_std_vector_x>::all_not<stan::is_std_vector>();
}
TEST(requires_prim_mat, any_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_std_vector_x>::any<stan::is_std_vector>();
}
TEST(requires_prim_mat, any_not_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_std_vector_x>::any_not<stan::is_std_vector>();
}

TEST(requires_prim_mat, eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_st, eigen_std_vector_x>::
      unary<std::is_floating_point>();
}
TEST(requires_prim_mat, not_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_st,
      eigen_std_vector_x>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_mat, all_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_st,
                            eigen_std_vector_x>::all<std::is_floating_point>();
}
TEST(requires_prim_mat, all_not_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_eigen_vector_st,
      eigen_std_vector_x>::all_not<std::is_floating_point>();
}
TEST(requires_prim_mat, any_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_st,
                            eigen_std_vector_x>::any<std::is_floating_point>();
}
TEST(requires_prim_mat, any_not_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_eigen_vector_st,
      eigen_std_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using std_vector_eigen_x = std::vector<Eigen::Matrix<T, -1, -1>>;

TEST(requires_prim_mat, std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt,
                            std_vector_eigen_x>::unary<stan::is_eigen>();
}
TEST(requires_prim_mat, not_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_vt,
                            std_vector_eigen_x>::not_unary<stan::is_eigen>();
}
TEST(requires_prim_mat, all_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt,
                            std_vector_eigen_x>::all<stan::is_eigen>();
}
TEST(requires_prim_mat, all_not_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt,
                            std_vector_eigen_x>::all_not<stan::is_eigen>();
}
TEST(requires_prim_mat, any_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt,
                            std_vector_eigen_x>::any<stan::is_eigen>();
}
TEST(requires_prim_mat, any_not_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt,
                            std_vector_eigen_x>::any_not<stan::is_eigen>();
}

TEST(requires_prim_mat, std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_st, std_vector_eigen_x>::
      unary<std::is_floating_point>();
}
TEST(requires_prim_mat, not_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_std_vector_st,
      std_vector_eigen_x>::not_unary<std::is_floating_point>();
}
TEST(requires_prim_mat, all_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_st,
                            std_vector_eigen_x>::all<std::is_floating_point>();
}
TEST(requires_prim_mat, all_not_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_std_vector_st,
      std_vector_eigen_x>::all_not<std::is_floating_point>();
}
TEST(requires_prim_mat, any_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_st,
                            std_vector_eigen_x>::any<std::is_floating_point>();
}
TEST(requires_prim_mat, any_not_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_std_vector_st,
      std_vector_eigen_x>::any_not<std::is_floating_point>();
}

TEST(requires_prim_mat, eigen_row_and_col) {
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::require_eigen_row_and_col_t;
  using stan::require_not_eigen_row_and_col_t;
  using stan::math::var;
  using stan::math::var_value;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_not_eigen_row_and_col_t,
                                        var_value<double>, VectorXd>::value));
  EXPECT_TRUE((require_variadic_checker<require_eigen_row_and_col_t,
                                        RowVectorXd, VectorXd>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_eigen_row_and_col_t,
                                         RowVectorXd, VectorXd>::value));

  EXPECT_FALSE((require_variadic_checker<require_eigen_row_and_col_t, VectorXd,
                                         VectorXd>::value));
  EXPECT_TRUE((require_variadic_checker<require_not_eigen_row_and_col_t,
                                        VectorXd, VectorXd>::value));
  EXPECT_FALSE((require_variadic_checker<require_eigen_row_and_col_t,
                                         RowVectorXd, RowVectorXd>::value));
  EXPECT_TRUE((require_variadic_checker<require_not_eigen_row_and_col_t,
                                        RowVectorXd, RowVectorXd>::value));
  EXPECT_FALSE((require_variadic_checker<require_eigen_row_and_col_t, VectorXd,
                                         RowVectorXd>::value));
  EXPECT_TRUE((require_variadic_checker<require_not_eigen_row_and_col_t,
                                        VectorXd, RowVectorXd>::value));
}
