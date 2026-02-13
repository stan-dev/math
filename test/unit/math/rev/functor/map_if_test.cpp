#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>
namespace {

TEST(MathPrim, map_if_base) {
  using stan::is_any_var_scalar;
  using stan::math::map_if;
  using stan::math::var;
  auto x
      = map_if<is_any_var_scalar>([](auto&& x) noexcept { return x + 1; },
                                        std::forward_as_tuple(1, var(2.0), std::vector<double>{3.0, 4.0}, std::vector<int>{4, 3, 2, 1}));
  static_assert(std::is_same_v<decltype(x), std::tuple<int, var, std::vector<double>, std::vector<int>>>,
                "Should be tuple of int, var, vector<double>, vector<int>!");
}

TEST(MathPrim, map_if_eigen) {
  using stan::math::map_if;
  using stan::test::is_ref_element_v;
  using stan::test::is_same_tuple_element_v;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  std::vector<double> b{1, 2, 3};
  auto x = map_if<stan::is_eigen>(
      [](auto&& x) { return (x * x).eval(); },
      std::forward_as_tuple(a, b, a * a, std::vector<int>{1, 2, 3}));
  EXPECT_MATRIX_EQ(std::get<0>(x), (a * a).eval());
  static_assert(is_same_tuple_element_v<0, decltype(x), Eigen::MatrixXd>,
                "1st should be MatrixXd!");
  static_assert(!is_ref_element_v<0, decltype(x)>, "0th should be an lvalue!");
  for (int i = 0; i < b.size(); i++) {
    EXPECT_EQ(std::get<1>(x)[i], b[i]);
  }
  static_assert(is_ref_element_v<1, decltype(x)>,
                "1st should be an reference!");
  static_assert(!is_ref_element_v<2, decltype(x)>, "2nd should be an lvalue!");
  static_assert(!is_ref_element_v<3, decltype(x)>, "3rd should be an lvalue!");
}
}  // namespace
