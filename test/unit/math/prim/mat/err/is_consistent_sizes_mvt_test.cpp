#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(isConsistentSizesMvt, isConsistentSizesMvt) {
  size_t expected_size = 3;
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

  EXPECT_TRUE(stan::math::is_consistent_sizes_mvt(good, good));
  EXPECT_TRUE(stan::math::is_consistent_sizes_mvt(good, good_single));
  EXPECT_TRUE(stan::math::is_consistent_sizes_mvt(good_single, good));
  EXPECT_TRUE(stan::math::is_consistent_sizes_mvt(good_single,
                                                      good_single));
  EXPECT_TRUE(stan::math::is_consistent_sizes_mvt(good_single,
                                                      bad_only_1));
  EXPECT_TRUE(stan::math::is_consistent_sizes_mvt(good_single,
                                                      bad_only_2));
  EXPECT_FALSE(stan::math::is_consistent_sizes_mvt(good, bad_only_2));
  EXPECT_FALSE(stan::math::is_consistent_sizes_mvt(good, bad_only_1));
  EXPECT_FALSE(stan::math::is_consistent_sizes_mvt(bad_only_2, good));
  EXPECT_FALSE(stan::math::is_consistent_sizes_mvt(bad_only_1, good));
}
