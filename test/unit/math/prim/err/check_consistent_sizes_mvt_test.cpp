#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(checkConsistentSizesMvt, checkConsistentSizesMvt) {
  using stan::math::check_consistent_sizes_mvt;

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
  std::vector<Eigen::VectorXd> empty;

  static const char* function
      = "TESTcheckConsistentSizesMvtcheckConsistentSizesMvt";

  EXPECT_NO_THROW(check_consistent_sizes_mvt(function));
  EXPECT_NO_THROW(check_consistent_sizes_mvt(function, "good", good));
  EXPECT_NO_THROW(
      check_consistent_sizes_mvt(function, "good_single", good_single));
  EXPECT_NO_THROW(check_consistent_sizes_mvt(function, "empty", empty));
  EXPECT_NO_THROW(check_consistent_sizes_mvt(function, "x1", good, "x2", good));
  EXPECT_NO_THROW(
      check_consistent_sizes_mvt(function, "x1", good, "x2", good_single));
  EXPECT_NO_THROW(
      check_consistent_sizes_mvt(function, "x1", good_single, "x2", good));
  EXPECT_NO_THROW(check_consistent_sizes_mvt(function, "x1", good_single, "x2",
                                             good_single));
  EXPECT_NO_THROW(check_consistent_sizes_mvt(function, "x1", good_single, "x2",
                                             bad_only_1));
  EXPECT_NO_THROW(check_consistent_sizes_mvt(function, "x1", good_single, "x2",
                                             bad_only_2));
  EXPECT_NO_THROW(check_consistent_sizes_mvt(function, "x1", good_single, "x2",
                                             bad_only_2, "x3", good_single,
                                             "x4", bad_only_2));
  EXPECT_THROW(
      check_consistent_sizes_mvt(function, "x1", good, "x2", bad_only_2),
      std::invalid_argument);
  EXPECT_THROW(
      check_consistent_sizes_mvt(function, "x1", good, "x2", bad_only_1),
      std::invalid_argument);
  EXPECT_THROW(
      check_consistent_sizes_mvt(function, "x1", bad_only_2, "x2", good),
      std::invalid_argument);
  EXPECT_THROW(
      check_consistent_sizes_mvt(function, "x1", bad_only_1, "x2", good),
      std::invalid_argument);
  EXPECT_THROW(
      check_consistent_sizes_mvt(function, "x1", good, "x2", good, "x3", good,
                                 "x4", good, "x5", bad_only_1, "x6", good),
      std::invalid_argument);
  EXPECT_THROW(check_consistent_sizes_mvt(function, "x1", good_single, "x2",
                                          good_single, "x3", good, "x4", good,
                                          "x5", bad_only_1, "x6", good),
               std::invalid_argument);
}
