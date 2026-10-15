#ifndef TEST_UNIT_MATH_MIX_FUN_MULTIPLY_OPERATOR_UTIL_HPP
#define TEST_UNIT_MATH_MIX_FUN_MULTIPLY_OPERATOR_UTIL_HPP

#include <test/unit/math/test_ad.hpp>

namespace stan {
namespace test {
namespace multiply_operator {

inline Eigen::MatrixXd test_matrix(int rows, int cols, double offset) {
  Eigen::MatrixXd x(rows, cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      x(i, j) = 0.3 * (i + 1) - 0.2 * (j + 1) + offset;
    }
  }
  return x;
}

}  // namespace multiply_operator
}  // namespace test
}  // namespace stan

#endif
