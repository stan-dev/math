#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat.hpp>
#include <test/unit/math/require_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <string>
#include <vector>

template <typename T>
using eigen_x = Eigen::Matrix<T, -1, -1>;

TEST(require_vt, eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vt,
                            eigen_x>::unary<std::is_floating_point>();
}
TEST(require_vt, not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vt,
                            eigen_x>::not_unary<std::is_floating_point>();
}
TEST(require_vt, all_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vt,
                            eigen_x>::all<std::is_floating_point>();
}
TEST(require_vt, all_not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vt,
                            eigen_x>::all_not<std::is_floating_point>();
}
TEST(require_vt, any_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vt,
                            eigen_x>::any<std::is_floating_point>();
}
TEST(require_vt, any_not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vt,
                            eigen_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_vector_x = Eigen::Matrix<T, 1, -1>;

TEST(require_vt, eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_vector_x>::unary<std::is_floating_point>();
}
TEST(require_vt, not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vector_vt, eigen_vector_x>::
      not_unary<std::is_floating_point>();
}
TEST(require_vt, all_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_vector_x>::all<std::is_floating_point>();
}
TEST(require_vt, all_not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_vector_x>::all_not<std::is_floating_point>();
}
TEST(require_vt, any_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_vector_x>::any<std::is_floating_point>();
}
TEST(require_vt, any_not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_eigen_x = Eigen::Matrix<Eigen::Matrix<T, -1, -1>, -1, -1>;

TEST(require_container_vt, eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vt,
                            eigen_eigen_x>::unary<stan::is_eigen>();
}
TEST(require_container_vt, not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vt,
                            eigen_eigen_x>::not_unary<stan::is_eigen>();
}
TEST(require_container_vt, all_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vt,
                            eigen_eigen_x>::all<stan::is_eigen>();
}
TEST(require_container_vt, all_not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vt,
                            eigen_eigen_x>::all_not<stan::is_eigen>();
}
TEST(require_container_vt, any_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vt,
                            eigen_eigen_x>::any<stan::is_eigen>();
}
TEST(require_container_vt, any_not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vt,
                            eigen_eigen_x>::any_not<stan::is_eigen>();
}

TEST(require_st, eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_st,
                            eigen_eigen_x>::unary<std::is_floating_point>();
}
TEST(require_st, not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_st,
                            eigen_eigen_x>::not_unary<std::is_floating_point>();
}
TEST(require_st, all_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_st,
                            eigen_eigen_x>::all<std::is_floating_point>();
}
TEST(require_st, all_not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_st,
                            eigen_eigen_x>::all_not<std::is_floating_point>();
}
TEST(require_st, any_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_st,
                            eigen_eigen_x>::any<std::is_floating_point>();
}
TEST(require_st, any_not_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_st,
                            eigen_eigen_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_eigen_vector_x = Eigen::Matrix<Eigen::Matrix<T, 1, -1>, 1, -1>;

TEST(require_container_vt, eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_eigen_vector_x>::unary<stan::is_vector>();
}
TEST(require_container_vt, not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::not_unary<stan::is_vector>();
}
TEST(require_container_vt, all_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_eigen_vector_x>::all<stan::is_vector>();
}
TEST(require_container_vt, all_not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::all_not<stan::is_vector>();
}
TEST(require_container_vt, any_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_eigen_vector_x>::any<stan::is_vector>();
}
TEST(require_container_vt, any_not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_eigen_vector_x>::any_not<stan::is_vector>();
}

TEST(require_st, eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_eigen_vector_st,
      eigen_eigen_vector_x>::unary<std::is_floating_point>();
}
TEST(require_st, not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_st,
      eigen_eigen_vector_x>::not_unary<std::is_floating_point>();
}
TEST(require_st, all_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_eigen_vector_st,
      eigen_eigen_vector_x>::all<std::is_floating_point>();
}
TEST(require_st, all_not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_eigen_vector_st,
      eigen_eigen_vector_x>::all_not<std::is_floating_point>();
}
TEST(require_st, any_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_eigen_vector_st,
      eigen_eigen_vector_x>::any<std::is_floating_point>();
}
TEST(require_st, any_not_eigen_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_eigen_vector_st,
      eigen_eigen_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using eigen_std_vector_x = Eigen::Matrix<std::vector<T>, 1, -1>;

TEST(require_container_vt, eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_vt,
                            eigen_std_vector_x>::unary<stan::is_std_vector>();
}
TEST(require_container_vt, not_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_vt,
      eigen_std_vector_x>::not_unary<stan::is_std_vector>();
}
TEST(require_container_vt, all_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_vt,
                            eigen_std_vector_x>::all<stan::is_std_vector>();
}
TEST(require_container_vt, all_not_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_eigen_vector_vt,
                            eigen_std_vector_x>::all_not<stan::is_std_vector>();
}
TEST(require_container_vt, any_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_vt,
                            eigen_std_vector_x>::any<stan::is_std_vector>();
}
TEST(require_container_vt, any_not_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_eigen_vector_vt,
                            eigen_std_vector_x>::any_not<stan::is_std_vector>();
}

TEST(require_st, eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_eigen_vector_st, eigen_std_vector_x>::
      unary<std::is_floating_point>();
}
TEST(require_st, not_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_eigen_vector_st,
      eigen_std_vector_x>::not_unary<std::is_floating_point>();
}
TEST(require_st, all_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_eigen_vector_st,
                            eigen_std_vector_x>::all<std::is_floating_point>();
}
TEST(require_st, all_not_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_eigen_vector_st,
      eigen_std_vector_x>::all_not<std::is_floating_point>();
}
TEST(require_st, any_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_eigen_vector_st,
                            eigen_std_vector_x>::any<std::is_floating_point>();
}
TEST(require_st, any_not_eigen_std_vector_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_eigen_vector_st,
      eigen_std_vector_x>::any_not<std::is_floating_point>();
}

template <typename T>
using std_vector_eigen_x = std::vector<Eigen::Matrix<T, -1, -1>>;

TEST(require_container_vt, std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_vt,
                            std_vector_eigen_x>::unary<stan::is_eigen>();
}
TEST(require_container_vt, not_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_not_std_vector_vt,
                            std_vector_eigen_x>::not_unary<stan::is_eigen>();
}
TEST(require_container_vt, all_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_vt,
                            std_vector_eigen_x>::all<stan::is_eigen>();
}
TEST(require_container_vt, all_not_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_not_std_vector_vt,
                            std_vector_eigen_x>::all_not<stan::is_eigen>();
}
TEST(require_container_vt, any_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_vt,
                            std_vector_eigen_x>::any<stan::is_eigen>();
}
TEST(require_container_vt, any_not_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_not_std_vector_vt,
                            std_vector_eigen_x>::any_not<stan::is_eigen>();
}

TEST(require_st, std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_std_vector_st, std_vector_eigen_x>::
      unary<std::is_floating_point>();
}
TEST(require_st, not_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_not_std_vector_st,
      std_vector_eigen_x>::not_unary<std::is_floating_point>();
}
TEST(require_st, all_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_all_std_vector_st,
                            std_vector_eigen_x>::all<std::is_floating_point>();
}
TEST(require_st, all_not_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_all_not_std_vector_st,
      std_vector_eigen_x>::all_not<std::is_floating_point>();
}
TEST(require_st, any_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<stan::require_any_std_vector_st,
                            std_vector_eigen_x>::any<std::is_floating_point>();
}
TEST(require_st, any_not_std_vector_eigen_test) {
  using stan::test::require_container_checker;
  require_container_checker<
      stan::require_any_not_std_vector_st,
      std_vector_eigen_x>::any_not<std::is_floating_point>();
}
