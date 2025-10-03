#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>
namespace {

TEST(MathMatrix, to_ref_matrix_exprs_tuple) {
  using stan::test::is_same_tuple_element_v;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  auto x = std::make_tuple(a * a, a, a.array() * 3);
  auto x_ref = stan::math::to_ref(x);
  using x_ref_t = decltype(x_ref);
  static_assert(is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>,
                "first entry should be Eigen::MatrixXd!");
  static_assert(is_same_tuple_element_v<1, x_ref_t, Eigen::MatrixXd>,
                "second entry should be Eigen::MatrixXd!");
  static_assert(
      is_same_tuple_element_v<2, x_ref_t,
                              Eigen::Array<double, -1, -1, 0, -1, -1>>,
      "third entry should be Eigen::ArrayXd!");
}

TEST(MathMatrix, to_ref_matrix_views_tuple) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  auto x = std::make_tuple(a.block(0, 0, 1, 1),
                           a(Eigen::all, std::vector{2, 1, 1}), a.array());
  auto x_ref = stan::math::to_ref(x);
  using x_ref_t = decltype(x_ref);
  using stan::test::is_same_tuple_element_v;
  static_assert(!is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>,
                "0th entry should be a view of an Eigen::MatrixXd!");
  static_assert(
      is_same_tuple_element_v<1, x_ref_t,
                              Eigen::Matrix<double, -1, -1, 0, -1, -1>>,
      "1st entry should be a view of an Eigen::MatrixXd!");
  static_assert(
      !is_same_tuple_element_v<2, x_ref_t,
                               Eigen::Array<double, -1, -1, 0, -1, -1>>,
      "2nd entry should be a view of an Eigen::ArrayXd!");
}

TEST(MathMatrix, to_ref_matrix_views_exprs_tuple) {
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  auto x = std::make_tuple(
      a.block(0, 0, 1, 1),
      std::make_tuple(a.block(0, 0, 1, 1), a(Eigen::all, std::vector{2, 1, 1}),
                      a.array()),
      std::make_tuple(a * a, a, a.array() * 3),
      a(Eigen::all, std::vector{2, 1, 1}), a.array() * a.array());
  auto x_ref = stan::math::to_ref(x);
  using x_ref_t = decltype(x_ref);
  using stan::test::is_same_tuple_element_v;
  static_assert(!is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>,
                "first entry should be a view of an Eigen::MatrixXd!");
  {
    using view_inner_tuple = std::tuple_element_t<1, decltype(x_ref)>;
    static_assert(
        !is_same_tuple_element_v<0, view_inner_tuple, Eigen::MatrixXd>,
        "tuple<1><0> entry should be a view of an Eigen::MatrixXd!");
    static_assert(
        is_same_tuple_element_v<1, view_inner_tuple,
                                Eigen::Matrix<double, -1, -1, 0, -1, -1>>,
        "tuple<1><1> entry should be Eigen::MatrixXd!");
    static_assert(
        !is_same_tuple_element_v<2, view_inner_tuple,
                                 Eigen::Array<double, -1, -1, 0, -1, -1>>,
        "tuple<1><2> entry should be a view of an Eigen::ArrayXd!");

    {
      using expr_inner_tuple = std::tuple_element_t<2, decltype(x_ref)>;
      static_assert(
          is_same_tuple_element_v<0, expr_inner_tuple, Eigen::MatrixXd>,
          "tuple<2><0> entry should be Eigen::MatrixXd!");
      static_assert(
          is_same_tuple_element_v<1, expr_inner_tuple, Eigen::MatrixXd>,
          "tuple<2><1> entry should be Eigen::MatrixXd!");
      static_assert(
          is_same_tuple_element_v<2, expr_inner_tuple,
                                  Eigen::Array<double, -1, -1, 0, -1, -1>>,
          "tuple<2><2> entry should be Eigen::ArrayXd!");
    }
    static_assert(
        is_same_tuple_element_v<3, x_ref_t,
                                Eigen::Matrix<double, -1, -1, 0, -1, -1>>,
        "tuple<3> entry should be Eigen::MatrixXd!");
    static_assert(
        is_same_tuple_element_v<4, x_ref_t,
                                Eigen::Array<double, -1, -1, 0, -1, -1>>,
        "tuple<4> entry should be Eigen::ArrayXd!");
  }
}
TEST(MathMatrix, to_ref_matrix_views_exprs_moves_tuple) {
  auto a = Eigen::MatrixXd::Random(3, 3);
  auto x_ref = stan::math::to_ref(std::forward_as_tuple(
      a.block(0, 0, 1, 1),
      std::forward_as_tuple(
          a.block(0, 0, 1, 1),
          a(Eigen::all,
            std::vector{Eigen::Index{2}, Eigen::Index{1}, Eigen::Index{1}}),
          a.array()),
      std::forward_as_tuple(a * a, a, a.array() * 3),
      a(Eigen::all,
        std::vector{Eigen::Index{2}, Eigen::Index{1}, Eigen::Index{1}}),
      a.array() * a.array()));
  using x_ref_t = decltype(x_ref);
  // These should all be evaluated
  using stan::test::is_same_tuple_element_v;
  static_assert(is_same_tuple_element_v<0, x_ref_t, Eigen::MatrixXd>,
                "first entry should be Eigen::MatrixXd!");
  {
    using view_inner_tuple = std::tuple_element_t<1, decltype(x_ref)>;
    static_assert(is_same_tuple_element_v<0, view_inner_tuple, Eigen::MatrixXd>,
                  "tuple<1><0> entry should be Eigen::MatrixXd!");
    static_assert(
        is_same_tuple_element_v<1, view_inner_tuple,
                                Eigen::Matrix<double, -1, -1, 0, -1, -1>>,
        "tuple<1><1> entry should be Eigen::MatrixXd!");
    static_assert(
        is_same_tuple_element_v<2, view_inner_tuple,
                                Eigen::Array<double, -1, -1, 0, -1, -1>>,
        "tuple<1><2> entry should be Eigen::ArrayXd!");
  }
  {
    using expr_inner_tuple = std::tuple_element_t<2, decltype(x_ref)>;
    static_assert(is_same_tuple_element_v<0, expr_inner_tuple, Eigen::MatrixXd>,
                  "tuple<2><0> entry should be Eigen::MatrixXd!");
    static_assert(is_same_tuple_element_v<1, expr_inner_tuple, Eigen::MatrixXd>,
                  "tuple<2><1> entry should be Eigen::MatrixXd!");
    static_assert(
        is_same_tuple_element_v<2, expr_inner_tuple,
                                Eigen::Array<double, -1, -1, 0, -1, -1>>,
        "tuple<2><2> entry should be Eigen::ArrayXd!");
  }
  static_assert(
      is_same_tuple_element_v<3, x_ref_t,
                              Eigen::Matrix<double, -1, -1, 0, -1, -1>>,
      "tuple<4> entry should be Eigen::ArrayXd!");
  static_assert(
      is_same_tuple_element_v<4, x_ref_t,
                              Eigen::Array<double, -1, -1, 0, -1, -1>>,
      "tuple<4> entry should be Eigen::ArrayXd!");
}
}  // namespace
