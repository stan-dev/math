#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

namespace holder_test {
template <typename T>
auto f(T&& a) {
  auto* a_heap = new std::remove_reference_t<T>(std::forward<T>(a));
  return stan::math::holder(*a_heap + *a_heap, a_heap);
}

template <typename T>
auto f2(T&& a) {
  auto* a_heap = new std::remove_reference_t<T>(std::forward<T>(a));
  return stan::math::holder(a_heap->array(), a_heap);
}

template <typename T>
auto f3(T&& a) {
  return stan::math::make_holder([](const auto& mat) { return mat + mat; },
                                 std::forward<T>(a));
}

template <typename T>
auto f4(T&& a) {
  return stan::math::make_holder([](auto& mat) { return mat.array(); },
                                 std::forward<T>(a));
}
}  // namespace holder_test

TEST(MathFunctions, holder_lvalue) {
  using Eigen::MatrixXd;
  MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd res = holder_test::f(m);
  MatrixXd correct = m + m;
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(MathFunctions, holder_rvalue) {
  using Eigen::MatrixXd;
  MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2 = m;
  MatrixXd res = holder_test::f(std::move(m2));
  MatrixXd correct = m + m;
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(MathFunctions, holder_lvalue_assign) {
  using Eigen::MatrixXd;
  MatrixXd m = MatrixXd::Identity(3, 3);
  MatrixXd m2(3, 3);
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  auto array_holder = holder_test::f2(m);
  array_holder = m2;
  MatrixXd res = array_holder;
  EXPECT_MATRIX_EQ(res, m2);
}

TEST(MathFunctions, make_holder_lvalue) {
  using Eigen::MatrixXd;
  MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  auto holder = holder_test::f3(m);
  m(2, 2) = 10;
  MatrixXd res = holder;
  MatrixXd correct = m + m;
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(MathFunctions, make_holder_rvalue) {
  using Eigen::MatrixXd;
  MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd m2 = m;
  MatrixXd res = holder_test::f3(std::move(m2));
  MatrixXd correct = m + m;
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(MathFunctions, make_holder_lvalue_assign) {
  using Eigen::MatrixXd;
  MatrixXd m = MatrixXd::Identity(3, 3);
  MatrixXd m2(3, 3);
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  auto array_holder = holder_test::f4(m);
  array_holder = m2;
  MatrixXd res = array_holder;
  EXPECT_MATRIX_EQ(res, m2);
  EXPECT_MATRIX_EQ(m, m2);
}

TEST(MathFunctions, make_holder_rvalue_assign) {
  using Eigen::MatrixXd;
  MatrixXd m = MatrixXd::Identity(3, 3);
  MatrixXd m2(3, 3);
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  auto array_holder = holder_test::f4(std::move(m));
  array_holder = m2;
  MatrixXd res = array_holder;
  EXPECT_MATRIX_EQ(res, m2);
}

TEST(MathFunctions, make_holder_assign_holder) {
  using Eigen::MatrixXd;
  MatrixXd m = MatrixXd::Identity(3, 3);
  MatrixXd m2(3, 3);
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  auto array_holder = holder_test::f4(m);
  array_holder = holder_test::f4(m2);
  MatrixXd res = array_holder;
  EXPECT_MATRIX_EQ(res, m2);
  EXPECT_MATRIX_EQ(m, m2);
}

TEST(MathFunctions, block_of_holder) {
  using Eigen::MatrixXd;
  MatrixXd m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  MatrixXd res = holder_test::f(m).block(1, 0, 2, 2);
  MatrixXd correct = (m + m).block(1, 0, 2, 2);
  EXPECT_MATRIX_EQ(res, correct);
}

TEST(MathFunctions, block_of_make_holder_assign) {
  using Eigen::MatrixXd;
  MatrixXd m = MatrixXd::Identity(3, 3);
  MatrixXd m2(3, 3);
  m2 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  auto array_holder = holder_test::f4(m);
  array_holder.block(0, 0, 3, 3) = m2;
  MatrixXd res = array_holder;
  EXPECT_MATRIX_EQ(res, m2);
  EXPECT_MATRIX_EQ(m, m2);
}
