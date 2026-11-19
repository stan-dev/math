#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>

// ref_type_t of an rvalue arena_matrix holds it by value, so a reference
// to a temporary that is destroyed cannot result
TEST_F(AgradRev, MathMetaRev_ref_type_rvalue_arena_matrix_by_value) {
  using stan::math::arena_matrix;
  using stan::math::var;
  EXPECT_FALSE((std::is_reference<
                stan::ref_type_t<arena_matrix<Eigen::VectorXd>&&>>::value));
  EXPECT_FALSE(
      (std::is_reference<
          stan::ref_type_t<arena_matrix<Eigen::Matrix<var, -1, 1>>&&>>::value));
  EXPECT_TRUE((std::is_reference<
               stan::ref_type_t<arena_matrix<Eigen::VectorXd>&>>::value));
}
