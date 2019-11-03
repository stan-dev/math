#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, is_mat_container_check) {
  using stan::is_container;
  EXPECT_TRUE(is_container<Eigen::MatrixXd>::value);
  EXPECT_TRUE(is_container<std::vector<double>>::value);
}
