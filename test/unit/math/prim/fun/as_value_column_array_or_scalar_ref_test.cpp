#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar_ref.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

TEST(MathFunPrim, as_value_column_array_or_scalar_ref_scalar) {
  using stan::math::as_value_column_array_or_scalar_ref;
  int a = 3;
  double b = 4;
  EXPECT_EQ(a, as_value_column_array_or_scalar_ref(a));
  EXPECT_EQ(b, as_value_column_array_or_scalar_ref(b));
}

TEST(MathFunPrim, as_value_column_array_or_scalar_ref_vectors) {
  using stan::math::as_value_column_array_or_scalar_ref;
  Eigen::VectorXd a(3);
  a << 1.5, -2, 3;
  Eigen::RowVectorXd b = a.transpose();
  std::vector<double> c{1.5, -2, 3};

  decltype(auto) a_val = as_value_column_array_or_scalar_ref(a);
  decltype(auto) b_val = as_value_column_array_or_scalar_ref(b);
  decltype(auto) c_val = as_value_column_array_or_scalar_ref(c);
  EXPECT_EQ(1, std::decay_t<decltype(b_val)>::ColsAtCompileTime);
  EXPECT_MATRIX_EQ(a.array(), a_val);
  EXPECT_MATRIX_EQ(a.array(), b_val);
  EXPECT_MATRIX_EQ(a.array(), c_val);
  // a plain vector is referenced, not copied
  a[1] = 7;
  EXPECT_EQ(7, a_val[1]);
}

TEST(MathFunPrim, as_value_column_array_or_scalar_ref_expression) {
  using stan::math::as_value_column_array_or_scalar_ref;
  Eigen::VectorXd a(3);
  a << 1.5, -2, 3;
  decltype(auto) a_val = as_value_column_array_or_scalar_ref(a * 2);
  // an expression is evaluated once
  EXPECT_TRUE(stan::is_plain_type<std::decay_t<decltype(a_val)>>::value);
  EXPECT_MATRIX_EQ((2 * a).array(), a_val);
}

TEST(MathFunPrim, as_value_column_array_or_scalar_ref_rvalue) {
  using stan::math::as_value_column_array_or_scalar_ref;
  Eigen::ArrayXd expected(3);
  expected << 1.5, -2, 3;
  // the result owns or copies the values of an rvalue argument
  decltype(auto) a_val
      = as_value_column_array_or_scalar_ref(std::vector<double>{1.5, -2, 3});
  decltype(auto) b_val
      = as_value_column_array_or_scalar_ref(Eigen::VectorXd(expected.matrix()));
  EXPECT_MATRIX_EQ(expected, a_val);
  EXPECT_MATRIX_EQ(expected, b_val);
}
