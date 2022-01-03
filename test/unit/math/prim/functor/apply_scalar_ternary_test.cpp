#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, apply_scalar_ternary) {
  Eigen::VectorXd vec1(2);
  Eigen::VectorXd res(2);
  vec1 << 1, 2;
  std::vector<double> stvec1{1, 2};

  // C, C, C
  Eigen::VectorXd y1 = stan::math::apply_scalar_ternary(
      vec1, vec1, vec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ(vec1 + vec1 + vec1, y1);
  std::vector<double> sty1 = stan::math::apply_scalar_ternary(
      stvec1, stvec1, stvec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ(vec1 + vec1 + vec1, stan::math::as_column_vector_or_scalar(sty1));

  // C, C, S
  y1 = stan::math::apply_scalar_ternary(
      vec1, vec1, 4,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + vec1.array() + 4).matrix(), y1);
  sty1 = stan::math::apply_scalar_ternary(
      stvec1, stvec1, 4,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + vec1.array() + 4).matrix(), stan::math::as_column_vector_or_scalar(sty1));

  // C, S, S
  y1 = stan::math::apply_scalar_ternary(
      vec1, 4, 4,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + 4 + 4).matrix(), y1);
  sty1 = stan::math::apply_scalar_ternary(
      stvec1, 4, 4,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + 4 + 4).matrix(), stan::math::as_column_vector_or_scalar(sty1));

  // C, S, C
  y1 = stan::math::apply_scalar_ternary(
      vec1, 4, vec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + 4 + vec1.array()).matrix(), y1);
  sty1 = stan::math::apply_scalar_ternary(
      stvec1, 4, stvec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + vec1.array() + 4).matrix(), stan::math::as_column_vector_or_scalar(sty1));

  // S, S, C
  y1 = stan::math::apply_scalar_ternary(
      4, 4, vec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((4 + 4 + vec1.array()).matrix(), y1);
  sty1 = stan::math::apply_scalar_ternary(
      4, 4, stvec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + 4 + 4).matrix(), stan::math::as_column_vector_or_scalar(sty1));

  // S, C, C
  y1 = stan::math::apply_scalar_ternary(
      4, vec1, vec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((4 + vec1.array() + vec1.array()).matrix(), y1);
  sty1 = stan::math::apply_scalar_ternary(
      4, stvec1, stvec1,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + vec1.array() + 4).matrix(), stan::math::as_column_vector_or_scalar(sty1));

  // S, C, S
  y1 = stan::math::apply_scalar_ternary(
      4, vec1, 4,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((4 + vec1.array() + 4).matrix(), y1);
  sty1 = stan::math::apply_scalar_ternary(
      4, stvec1, 4,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_MATRIX_FLOAT_EQ((vec1.array() + 4 + 4).matrix(), stan::math::as_column_vector_or_scalar(sty1));

  // S, S, S
  auto y2 = stan::math::apply_scalar_ternary(
      4, 4, 4,
      [&](const auto& x, const auto& y, const auto& z) { return x + y + z; });
  EXPECT_FLOAT_EQ((4 + 4 + 4), y2);
}
