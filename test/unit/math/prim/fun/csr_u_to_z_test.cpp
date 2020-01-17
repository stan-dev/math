#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <exception>

TEST(MathUtoZ, sizeZeroInput) {
  using stan::math::csr_u_to_z;
  std::vector<int> u(0);
  for (int n = -1; n < 2; ++n)
    EXPECT_THROW(csr_u_to_z(u, n), std::out_of_range);
}
TEST(MathUtoZ, nonemptyInput) {
  using stan::math::csr_u_to_z;
  std::vector<int> u(3);
  EXPECT_THROW(csr_u_to_z(u, 10), std::out_of_range);
  EXPECT_NO_THROW(csr_u_to_z(u, 0));
  EXPECT_THROW(csr_u_to_z(u, -2), std::out_of_range);
}
TEST(MathUtoZ, outputs) {
  using stan::math::csr_u_to_z;
  std::vector<int> u;
  u.push_back(3);
  u.push_back(4);
  u.push_back(9);
  EXPECT_EQ(1, csr_u_to_z(u, 0));
  EXPECT_EQ(5, csr_u_to_z(u, 1));
}
