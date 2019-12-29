#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, is_scalar_container_check) {
  using stan::is_container;
  EXPECT_FALSE(is_container<int>::value);
  EXPECT_FALSE(is_container<double>::value);
}

TEST(MetaTraits, is_mat_container_check) {
  using stan::is_container;
  EXPECT_TRUE(is_container<Eigen::MatrixXd>::value);
  EXPECT_TRUE(is_container<std::vector<double>>::value);
}

TEST(MetaTraits, is_container_container_check) {
  using stan::is_container;
  EXPECT_TRUE(is_container<std::vector<Eigen::MatrixXd>>::value);
  EXPECT_TRUE(is_container<std::vector<std::vector<double>>>::value);
}
