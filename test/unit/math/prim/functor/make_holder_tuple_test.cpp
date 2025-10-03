#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>
namespace {

TEST(MathMatrix, partial_forward_as_tuples_rvalues) {
  using stan::test::is_const_ref_element_v;
  using stan::test::is_lvalue_ref_element_v;
  using stan::test::is_ref_element_v;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  const Eigen::MatrixXd a_const = Eigen::MatrixXd::Random(3, 3);
  auto a_expr = Eigen::MatrixXd::Random(3, 3);
  auto x_fwd_tuple = stan::math::make_holder_tuple(
      a * a, a, a.array() * 3, Eigen::MatrixXd::Random(3, 3), a_const,
      a_const * a_const, a_expr, a_expr * a_expr);
  using x_ref_t = decltype(x_fwd_tuple);
  static_assert(!is_ref_element_v<0, x_ref_t>,
                "0th entry should be an lvalue!");
  static_assert(is_lvalue_ref_element_v<1, x_ref_t>,
                "1st entry should be an lvalue reference!");
  static_assert(!is_ref_element_v<2, x_ref_t>,
                "2nd entry should be an lvalue!");
  static_assert(!is_ref_element_v<3, x_ref_t>,
                "3rd entry should be an lvalue!");
  static_assert(is_const_ref_element_v<4, x_ref_t>,
                "4th entry should be an const lvalue reference!");
  static_assert(!is_ref_element_v<5, x_ref_t>,
                "5th entry should be an lvalue!");
  static_assert(is_lvalue_ref_element_v<6, x_ref_t>,
                "6th entry should be an lvalue reference!");
  static_assert(!is_ref_element_v<7, x_ref_t>,
                "7th entry should be an lvalue!");
}
}  // namespace
