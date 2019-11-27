#ifndef TEST_UNIT_MATH_TEST_IS_FINITE_HPP
#define TEST_UNIT_MATH_TEST_IS_FINITE_HPP

#include <stan/math.hpp>
#include <vector>

namespace stan {
namespace test {

/**
 * Return true if the specified value is finite.
 *
 * @param x value to test
 * @return true if value is finite
 */
bool is_finite(double x) {
  return !stan::math::is_inf(x) && !stan::math::is_nan(x);
}

/**
 * Return true if all of the elements in the container are finite
 *
 * @tparam T scalar type
 * @tparam R row type
 * @tparam C col type
 * @param x container to test
 * @return true if all container values are finite
 */
template <typename T, int R, int C>
bool is_finite(const Eigen::Matrix<T, R, C>& x) {
  for (int i = 0; i < x.size(); ++i)
    if (!is_finite(x(i)))
      return false;
  return true;
}

/**
 * Return true if all of the elements in the container are finite
 *
 * @tparam T contained type
 * @param x container to test
 * @return true if all container values are finite
 */
template <typename T>
bool is_finite(const std::vector<T>& x) {
  for (size_t i = 0; i < x.size(); ++i)
    if (!is_finite(x[i]))
      return false;
  return true;
}

}  // namespace test
}  // namespace stan
#endif
