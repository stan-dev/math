#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>
#include <vector>

template <typename T>
using eigen_x = Eigen::Matrix<T, -1, -1>;

// Test same
TEST(requires, same_vtest) {
  using stan::require_same_vt;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<stan::require_same_vt, eigen_x<double>,
                                         eigen_x<int>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_same_vt, eigen_x<double>,
                                        eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_same_vt, eigen_x<int>,
                                        eigen_x<int>>::value));
}

TEST(requires, not_same_vtest) {
  using stan::require_not_same_vt;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_not_same_vt, eigen_x<double>,
                                        eigen_x<int>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_same_vt, eigen_x<double>,
                                         eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_same_vt, eigen_x<int>,
                                         eigen_x<int>>::value));
}

TEST(requires, all_same_vtest) {
  using stan::require_all_same_vt;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_same_vt, eigen_x<double>,
                                         std::string, eigen_x<double>>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_all_same_vt, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_same_vt, eigen_x<int>,
                                        eigen_x<int>, eigen_x<int>>::value));
}

TEST(requires, all_not_same_vtest) {
  using stan::require_any_not_same_vt;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE(
      (require_variadic_checker<require_any_not_same_vt, eigen_x<double>,
                                eigen_x<int>, eigen_x<double>>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_any_not_same_vt, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_same_vt, eigen_x<int>,
                                         eigen_x<int>, eigen_x<int>>::value));
}

// Test same
TEST(requires, same_stest) {
  using stan::require_same_st;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<stan::require_same_st, eigen_x<double>,
                                         eigen_x<int>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_same_st, eigen_x<double>,
                                        eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<stan::require_same_st, eigen_x<int>,
                                        eigen_x<int>>::value));
}

TEST(requires, not_same_stest) {
  using stan::require_not_same_st;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE((require_variadic_checker<require_not_same_st, eigen_x<double>,
                                        eigen_x<int>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_same_st, eigen_x<double>,
                                         eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_not_same_st, eigen_x<int>,
                                         eigen_x<int>>::value));
}

TEST(requires, all_same_stest) {
  using stan::require_all_same_st;
  using stan::test::require_variadic_checker;
  EXPECT_FALSE((require_variadic_checker<require_all_same_st, eigen_x<double>,
                                         std::string, eigen_x<double>>::value));
  EXPECT_TRUE(
      (require_variadic_checker<require_all_same_st, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_TRUE((require_variadic_checker<require_all_same_st, eigen_x<int>,
                                        eigen_x<int>, eigen_x<int>>::value));
}

TEST(requires, all_not_same_stest) {
  using stan::require_any_not_same_st;
  using stan::test::require_variadic_checker;
  EXPECT_TRUE(
      (require_variadic_checker<require_any_not_same_st, eigen_x<double>,
                                eigen_x<int>, eigen_x<double>>::value));
  EXPECT_FALSE(
      (require_variadic_checker<require_any_not_same_st, eigen_x<double>,
                                eigen_x<double>, eigen_x<double>>::value));
  EXPECT_FALSE((require_variadic_checker<require_any_not_same_st, eigen_x<int>,
                                         eigen_x<int>, eigen_x<int>>::value));
}

TEST(requires, eigen_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_eigen_t, eigen_x<double>>::unary();
}
TEST(requires, eigen_not_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_not_eigen_t, eigen_x<double>,
                       eigen_x<double>>::not_unary();
}
TEST(requires, eigen_all_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_eigen_t, eigen_x<double>,
                       eigen_x<double>>::all();
}
TEST(requires, eigen_all_not_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_all_not_eigen_t, eigen_x<double>,
                       eigen_x<double>>::all_not();
}
TEST(requires, eigen_any_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_eigen_t, eigen_x<double>,
                       eigen_x<double>>::any();
}
TEST(requires, eigen_any_not_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_any_not_eigen_t, eigen_x<double>,
                       eigen_x<double>>::any_not();
}

TEST(requires, eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vt,
                            eigen_x>::unary<std::is_floating_point>();
}
TEST(requires, not_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vt,
                            eigen_x>::not_unary<std::is_floating_point>();
}
TEST(requires, all_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vt,
                            eigen_x>::all<std::is_floating_point>();
}
TEST(requires, all_not_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vt,
                            eigen_x>::all_not<std::is_floating_point>();
}
TEST(requires, any_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vt,
                            eigen_x>::any<std::is_floating_point>();
}
TEST(requires, any_not_eigen_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vt,
                            eigen_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_vector_x = Eigen::Matrix<T, 1, -1>;

TEST(requires, eigen_vector_t_test) {
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

TEST(requires, eigen_std_vector_t_test) {
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

TEST(requires, eigen_vector_like_t_test) {
  using stan::test::require_scal_checker;
  require_scal_checker<stan::require_vector_like_t, eigen_x<double>,
                       eigen_vector_x<double>>::unary();
  require_scal_checker<stan::require_not_vector_like_t, eigen_x<double>,
                       eigen_vector_x<double>>::not_unary();
  require_scal_checker<stan::require_all_vector_like_t, eigen_x<double>,
                       eigen_vector_x<double>>::all();
  require_scal_checker<stan::require_all_not_vector_like_t, eigen_x<double>,
                       eigen_vector_x<double>>::all_not();
  require_scal_checker<stan::require_any_vector_like_t, eigen_x<double>,
                       eigen_vector_x<double>>::any();
  require_scal_checker<stan::require_any_not_vector_like_t, std::vector<double>,
                       std::vector<double>>::any_not();
}

TEST(requires, eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_vector_x>::unary<std::is_floating_point>();
}
TEST(requires, not_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vector_vt, eigen_vector_x>::
      not_unary<std::is_floating_point>();
}
TEST(requires, all_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_vector_x>::all<std::is_floating_point>();
}
TEST(requires, all_not_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_vector_x>::all_not<std::is_floating_point>();
}
TEST(requires, any_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_vector_x>::any<std::is_floating_point>();
}
TEST(requires, any_not_eigen_vector_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_eigen_x = Eigen::Matrix<Eigen::Matrix<T, -1, -1>, -1, -1>;

TEST(requires, eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vt,
                            eigen_eigen_x>::unary<stan::is_eigen>();
}
TEST(requires, not_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vt,
                            eigen_eigen_x>::not_unary<stan::is_eigen>();
}
TEST(requires, all_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vt,
                            eigen_eigen_x>::all<stan::is_eigen>();
}
TEST(requires, all_not_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vt,
                            eigen_eigen_x>::all_not<stan::is_eigen>();
}
TEST(requires, any_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vt,
                            eigen_eigen_x>::any<stan::is_eigen>();
}
TEST(requires, any_not_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vt,
                            eigen_eigen_x>::any_not<stan::is_eigen>();
}

TEST(requires, eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_st,
                            eigen_eigen_x>::unary<std::is_floating_point>();
}
TEST(requires, not_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_st,
                            eigen_eigen_x>::not_unary<std::is_floating_point>();
}
TEST(requires, all_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_st,
                            eigen_eigen_x>::all<std::is_floating_point>();
}
TEST(requires, all_not_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_st,
                            eigen_eigen_x>::all_not<std::is_floating_point>();
}
TEST(requires, any_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_st,
                            eigen_eigen_x>::any<std::is_floating_point>();
}
TEST(requires, any_not_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_st,
                            eigen_eigen_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_eigen_vector_x = Eigen::Matrix<Eigen::Matrix<T, 1, -1>, 1, -1>;

TEST(requires, eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_eigen_vector_x>::unary<stan::is_vector>();
}
TEST(requires, not_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::not_unary<stan::is_vector>();
}
TEST(requires, all_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_eigen_vector_x>::all<stan::is_vector>();
}
TEST(requires, all_not_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::all_not<stan::is_vector>();
}
TEST(requires, any_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_eigen_vector_x>::any<stan::is_vector>();
}
TEST(requires, any_not_eigen_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::any_not<stan::is_vector>();
}

TEST(requires, eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_eigen_vector_st,
      eigen_eigen_vector_x>::unary<std::is_floating_point>();
}
TEST(requires, not_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_st,
      eigen_eigen_vector_x>::not_unary<std::is_floating_point>();
}
TEST(requires, all_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_eigen_vector_st,
      eigen_eigen_vector_x>::all<std::is_floating_point>();
}
TEST(requires, all_not_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_eigen_vector_st,
      eigen_eigen_vector_x>::all_not<std::is_floating_point>();
}
TEST(requires, any_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_eigen_vector_st,
      eigen_eigen_vector_x>::any<std::is_floating_point>();
}
TEST(requires, any_not_eigen_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_eigen_vector_st,
      eigen_eigen_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_std_vector_x = Eigen::Matrix<std::vector<T>, 1, -1>;

TEST(requires, eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_std_vector_x>::unary<stan::is_std_vector>();
}
TEST(requires, not_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_vt,
      eigen_std_vector_x>::not_unary<stan::is_std_vector>();
}
TEST(requires, all_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_std_vector_x>::all<stan::is_std_vector>();
}
TEST(requires, all_not_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_std_vector_x>::all_not<stan::is_std_vector>();
}
TEST(requires, any_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_std_vector_x>::any<stan::is_std_vector>();
}
TEST(requires, any_not_eigen_std_vector_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_std_vector_x>::any_not<stan::is_std_vector>();
}

TEST(requires, eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_st, eigen_std_vector_x>::
      unary<std::is_floating_point>();
}
TEST(requires, not_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_st,
      eigen_std_vector_x>::not_unary<std::is_floating_point>();
}
TEST(requires, all_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_st,
                            eigen_std_vector_x>::all<std::is_floating_point>();
}
TEST(requires, all_not_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_eigen_vector_st,
      eigen_std_vector_x>::all_not<std::is_floating_point>();
}
TEST(requires, any_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_st,
                            eigen_std_vector_x>::any<std::is_floating_point>();
}
TEST(requires, any_not_eigen_std_vector_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_eigen_vector_st,
      eigen_std_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using std_vector_eigen_x = std::vector<Eigen::Matrix<T, -1, -1>>;

TEST(requires, std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt,
                            std_vector_eigen_x>::unary<stan::is_eigen>();
}
TEST(requires, not_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_vt,
                            std_vector_eigen_x>::not_unary<stan::is_eigen>();
}
TEST(requires, all_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt,
                            std_vector_eigen_x>::all<stan::is_eigen>();
}
TEST(requires, all_not_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt,
                            std_vector_eigen_x>::all_not<stan::is_eigen>();
}
TEST(requires, any_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt,
                            std_vector_eigen_x>::any<stan::is_eigen>();
}
TEST(requires, any_not_std_vector_eigen_container_vt_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt,
                            std_vector_eigen_x>::any_not<stan::is_eigen>();
}

TEST(requires, std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_st, std_vector_eigen_x>::
      unary<std::is_floating_point>();
}
TEST(requires, not_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_std_vector_st,
      std_vector_eigen_x>::not_unary<std::is_floating_point>();
}
TEST(requires, all_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_st,
                            std_vector_eigen_x>::all<std::is_floating_point>();
}
TEST(requires, all_not_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_std_vector_st,
      std_vector_eigen_x>::all_not<std::is_floating_point>();
}
TEST(requires, any_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_st,
                            std_vector_eigen_x>::any<std::is_floating_point>();
}
TEST(requires, any_not_std_vector_eigen_st_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_std_vector_st,
      std_vector_eigen_x>::any_not<std::is_floating_point>();
}
