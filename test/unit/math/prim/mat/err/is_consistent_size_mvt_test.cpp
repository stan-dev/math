#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(isConsistentSizeMvt, isConsistentSizeMvt) {
  size_t size = 3;

  std::vector<Eigen::VectorXd> good(3);
  Eigen::VectorXd good_single(4);
  good_single << 1.0, 2.0, 3.0, 4.0;
  good[0] = good_single;
  good[1] = good_single;
  good[2] = good_single;
  std::vector<Eigen::VectorXd> bad_only_2(2);
  bad_only_2[0] = good_single;
  bad_only_2[1] = good_single;
  std::vector<Eigen::VectorXd> bad_only_1(1);
  bad_only_1[0] = good_single;
  std::vector<Eigen::VectorXd> empty;

  EXPECT_TRUE(stan::math::is_consistent_size_mvt(good, size));
  EXPECT_TRUE(stan::math::is_consistent_size_mvt(good_single, size));
  EXPECT_TRUE(stan::math::is_consistent_size_mvt(empty, 0));
  EXPECT_FALSE(stan::math::is_consistent_size_mvt(bad_only_2, size));
  EXPECT_FALSE(stan::math::is_consistent_size_mvt(bad_only_1, size));
  EXPECT_FALSE(stan::math::is_consistent_size_mvt(empty, size));
}
