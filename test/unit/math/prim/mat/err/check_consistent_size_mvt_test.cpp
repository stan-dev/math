#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(checkConsistentSizeMvt, checkConsistentSizeMvt) {
  using stan::math::check_consistent_size_mvt;

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
      = "TESTcheckConsistentSizeMvtcheckConsistentSizeMvt";

  EXPECT_NO_THROW(
      check_consistent_size_mvt(function, "good", good, expected_size));
  EXPECT_NO_THROW(check_consistent_size_mvt(function, "good_single",
                                            good_single, expected_size));
  EXPECT_NO_THROW(check_consistent_size_mvt(function, "empty", empty, 0));
  EXPECT_THROW(check_consistent_size_mvt(function, "bad_only_2", bad_only_2,
                                         expected_size),
               std::invalid_argument);
  EXPECT_THROW(check_consistent_size_mvt(function, "bad_only_1", bad_only_1,
                                         expected_size),
               std::invalid_argument);
  EXPECT_THROW(
      check_consistent_size_mvt(function, "empty", empty, expected_size),
      std::invalid_argument);
}
