#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>
namespace {

template <Eigen::Index Idx, typename Tuple, typename T>
static constexpr bool is_same_tuple_element_v = std::is_same<std::decay_t<std::tuple_element_t<Idx, Tuple>>, T>::value;

TEST(MathMatrix, to_ref_matrix_exprs_tuple) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  auto x = std::make_tuple(a * a, a, a.array() * 3);
  auto x_ref = stan::math::to_ref(x);
  std::cout << stan::math::test::type_name<decltype(x)>() << std::endl;
  std::cout << stan::math::test::type_name<decltype(x_ref)>() << std::endl;
  using x_ref_t = decltype(x_ref);
  static_assert(is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>, "first entry should be Eigen::MatrixXd!");
  static_assert(is_same_tuple_element_v<1, x_ref_t, Eigen::MatrixXd>, "second entry should be Eigen::MatrixXd!");
  static_assert(is_same_tuple_element_v<2, x_ref_t, Eigen::Array<double, -1, -1, 0, -1, -1>>, "third entry should be Eigen::ArrayXd!");
}

TEST(MathMatrix, to_ref_matrix_views_tuple) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  auto x = std::make_tuple(a.block(0, 0, 1, 1), a(Eigen::placeholders::all,{2, 1, 1}), a.array());
  auto x_ref = stan::math::to_ref(x);
  std::cout << stan::math::test::type_name<decltype(x)>() << std::endl;
  std::cout << stan::math::test::type_name<decltype(x_ref)>() << std::endl;
  using x_ref_t = decltype(x_ref);
  static_assert(!is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>, "first entry should be Eigen::MatrixXd!");
  static_assert(!is_same_tuple_element_v<1, x_ref_t, Eigen::MatrixXd>, "second entry should be Eigen::MatrixXd!");
  static_assert(!is_same_tuple_element_v<2, x_ref_t, Eigen::Array<double, -1, -1, 0, -1, -1>>, "third entry should be Eigen::MatrixXd!");
}

TEST(MathMatrix, to_ref_matrix_views_exprs_tuple) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  auto x = std::make_tuple(
    a.block(0, 0, 1, 1),
    std::make_tuple(
      a.block(0, 0, 1, 1),
      a(Eigen::placeholders::all,{2, 1, 1}),
      a.array()),
    std::make_tuple(
      a * a,
      a,
      a.array() * 3),
    a(Eigen::placeholders::all,{2, 1, 1}),
    a.array() * a.array());
  auto x_ref = stan::math::to_ref(x);
  using x_ref_t = decltype(x_ref);
  static_assert(!is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>, "first entry should be Eigen::MatrixXd!");
  {
    using view_inner_tuple = std::tuple_element_t<1, decltype(x_ref)>;
    static_assert(!is_same_tuple_element_v<0, view_inner_tuple, Eigen::MatrixXd>, "tuple<1><0> entry should be Eigen::MatrixXd!");
    // TODO(Steve): Views like this should still be unevaluated
    static_assert(is_same_tuple_element_v<1, view_inner_tuple, Eigen::Matrix<double, -1, 3, 0, -1, 3>>, "tuple<1><1> entry should be Eigen::MatrixXd!");
    static_assert(!is_same_tuple_element_v<2, view_inner_tuple, Eigen::Array<double, -1, -1, 0, -1, -1>>, "tuple<1><2> entry should be Eigen::MatrixXd!");

  {
    using expr_inner_tuple = std::tuple_element_t<2, decltype(x_ref)>;
    static_assert(is_same_tuple_element_v<0, expr_inner_tuple, Eigen::MatrixXd>, "tuple<2><0> entry should be Eigen::MatrixXd!");
    static_assert(is_same_tuple_element_v<1, expr_inner_tuple, Eigen::MatrixXd>, "tuple<2><1> entry should be Eigen::MatrixXd!");
    static_assert(is_same_tuple_element_v<2, expr_inner_tuple, Eigen::Array<double, -1, -1, 0, -1, -1>>, "tuple<2><2> entry should be Eigen::MatrixXd!");
  }

  static_assert(is_same_tuple_element_v<3, x_ref_t, Eigen::Matrix<double, -1, 3, 0, -1, 3>>, "tuple<4> entry should be Eigen::ArrayXd!");
  static_assert(is_same_tuple_element_v<4, x_ref_t, Eigen::Array<double, -1, -1, 0, -1, -1>>, "tuple<4> entry should be Eigen::ArrayXd!");
}
}
TEST(MathMatrix, to_ref_matrix_views_exprs_moves_tuple) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  auto x = std::tuple(
    Eigen::MatrixXd::Random(3, 3).block(0, 0, 1, 1),
    std::forward_as_tuple(
      Eigen::MatrixXd::Random(3, 3).block(0, 0, 1, 1),
      Eigen::MatrixXd::Random(3, 3)(Eigen::placeholders::all,{Eigen::Index{2}, Eigen::Index{1}, Eigen::Index{1}}),
      Eigen::MatrixXd::Random(3, 3).array()),
    std::forward_as_tuple(
      Eigen::MatrixXd::Random(3, 3) * Eigen::MatrixXd::Random(3, 3),
      Eigen::MatrixXd::Random(3, 3),
      Eigen::MatrixXd::Random(3, 3).array() * 3),
    Eigen::MatrixXd::Random(3, 3)(Eigen::placeholders::all,{Eigen::Index{2}, Eigen::Index{1}, Eigen::Index{1}}),
    Eigen::MatrixXd::Random(3, 3).array() * Eigen::MatrixXd::Random(3, 3).array());
  std::cout << "orig: \n";
  std::cout << stan::math::test::type_name<decltype(x)>() << std::endl;
  using x_ref_t = decltype(stan::math::to_ref(std::move(x)));
  std::cout << "move: \n";
  std::cout << stan::math::test::type_name<x_ref_t>() << std::endl;
  auto x_ref = stan::math::to_ref(std::move(x));
  /*
  static_assert(!is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>, "first entry should be Eigen::MatrixXd!");
  {
    using view_inner_tuple = std::tuple_element_t<1, decltype(x_ref)>;
    static_assert(!is_same_tuple_element_v<0, view_inner_tuple, Eigen::MatrixXd>, "tuple<1><0> entry should be Eigen::MatrixXd!");
    // TODO(Steve): Views like this should still be unevaluated
    static_assert(is_same_tuple_element_v<1, view_inner_tuple, Eigen::Matrix<double, -1, 3, 0, -1, 3>>, "tuple<1><1> entry should be Eigen::MatrixXd!");
    static_assert(!is_same_tuple_element_v<2, view_inner_tuple, Eigen::Array<double, -1, -1, 0, -1, -1>>, "tuple<1><2> entry should be Eigen::MatrixXd!");

  {
    using expr_inner_tuple = std::tuple_element_t<2, decltype(x_ref)>;
    static_assert(is_same_tuple_element_v<0, expr_inner_tuple, Eigen::MatrixXd>, "tuple<2><0> entry should be Eigen::MatrixXd!");
    static_assert(is_same_tuple_element_v<1, expr_inner_tuple, Eigen::MatrixXd>, "tuple<2><1> entry should be Eigen::MatrixXd!");
    static_assert(is_same_tuple_element_v<2, expr_inner_tuple, Eigen::Array<double, -1, -1, 0, -1, -1>>, "tuple<2><2> entry should be Eigen::MatrixXd!");
  }

  static_assert(is_same_tuple_element_v<3, x_ref_t, Eigen::Matrix<double, -1, 3, 0, -1, 3>>, "tuple<4> entry should be Eigen::ArrayXd!");
  static_assert(is_same_tuple_element_v<4, x_ref_t, Eigen::Array<double, -1, -1, 0, -1, -1>>, "tuple<4> entry should be Eigen::ArrayXd!");
  */
}
}


