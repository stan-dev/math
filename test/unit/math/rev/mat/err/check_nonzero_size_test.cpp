#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <limits>
#include <vector>

TEST(AgradRevErrorHandlingMatrix, checkNonzeroSizeMatrix) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> y;
  using stan::math::check_nonzero_size;
  var result;

  y.resize(3, 3);
  EXPECT_NO_THROW(check_nonzero_size("checkNonzeroSize", "y", y));
  y.resize(2, 3);
  EXPECT_NO_THROW(check_nonzero_size("checkNonzeroSize", "y", y));

  y.resize(0, 0);
  EXPECT_THROW_MSG(check_nonzero_size("checkNonzeroSize", "y", y),
                   std::invalid_argument, "y has size 0");

  std::vector<var> a;
  a.push_back(3.0);
  a.push_back(3.0);
  a.push_back(3.0);
  a.push_back(3.0);


  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize",
                                                 "a", a));

  a.resize(2);
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize",
                                                 "a", a));

  a.resize(0);
  EXPECT_THROW_MSG(stan::math::check_nonzero_size("checkNonzeroSize", "a", a),
                   std::invalid_argument,
                   "a has size 0");
}

TEST(AgradRevErrorHandlingMatrix, checkNonzeroSizeMatrix_nan) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> y;
  var result;
  var nan = std::numeric_limits<var>::quiet_NaN();

  y.resize(3, 3);
  y << nan, nan, nan, nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize",
                                                 "y", y));
  y.resize(2, 3);
  y << nan, nan, nan, nan, nan, nan;
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize",
                                                 "y", y));

  y.resize(0, 0);
  EXPECT_THROW_MSG(stan::math::check_nonzero_size("checkNonzeroSize", "y", y),
                   std::invalid_argument,
                   "has size 0");

  std::vector<var> a;
  a.push_back(nan);
  a.push_back(nan);
  a.push_back(nan);
  a.push_back(nan);

  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize",
                                                 "a", a));

  a.resize(2);
  EXPECT_NO_THROW(stan::math::check_nonzero_size("checkNonzeroSize",
                                                 "a", a));

  a.resize(0);
  EXPECT_THROW_MSG(stan::math::check_nonzero_size("checkNonzeroSize", "a", a),
                   std::invalid_argument,
                   "a has size 0");
}
