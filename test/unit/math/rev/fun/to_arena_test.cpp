#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRev, to_arena_scalar_test) {
  int a = 2;
  auto b = stan::math::to_arena(a);
  EXPECT_EQ(b, a);
  EXPECT_TRUE((std::is_same<decltype(a), decltype(b)>::value));
  double a2 = 2;
  auto b2 = stan::math::to_arena(a2);
  EXPECT_EQ(b2, a2);
  EXPECT_TRUE((std::is_same<decltype(a2), decltype(b2)>::value));
}

TEST(AgradRev, to_arena_std_vector_test) {
  std::vector<int> a{1, 2};
  auto b = stan::math::to_arena(a);
  ASSERT_EQ(a.size(), b.size());
  for (int i = 0; i < a.size(); i++) {
    EXPECT_EQ(a[i], b[i]);
  }
  EXPECT_FALSE((std::is_same<decltype(a), decltype(b)>::value));

  auto c = stan::math::to_arena(b);
  EXPECT_EQ(b.size(), c.size());
  EXPECT_EQ(b.data(), c.data());
}

TEST(AgradRev, to_arena_col_vector_test) {
  Eigen::VectorXd a(2);
  a << 1, 2;
  auto b = stan::math::to_arena(a);
  EXPECT_MATRIX_EQ(a, b);
  EXPECT_FALSE((std::is_same<decltype(a), decltype(b)>::value));
  auto c = stan::math::to_arena(b);
  EXPECT_EQ(b.size(), c.size());
  EXPECT_EQ(b.data(), c.data());
}

TEST(AgradRev, to_arena_row_vector_test) {
  Eigen::RowVectorXd a(2);
  a << 1, 2;
  auto b = stan::math::to_arena(a);
  EXPECT_MATRIX_EQ(a, b);
  EXPECT_FALSE((std::is_same<decltype(a), decltype(b)>::value));
  auto c = stan::math::to_arena(b);
  EXPECT_EQ(b.size(), c.size());
  EXPECT_EQ(b.data(), c.data());
}

TEST(AgradRev, to_arena_matrix_test) {
  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  auto b = stan::math::to_arena(a);
  EXPECT_MATRIX_EQ(a, b);
  EXPECT_FALSE((std::is_same<decltype(a), decltype(b)>::value));
  auto c = stan::math::to_arena(b);
  EXPECT_EQ(b.size(), c.size());
  EXPECT_EQ(b.data(), c.data());
}
