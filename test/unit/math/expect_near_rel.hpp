#ifndef TEST_UNIT_MATH_TEST_EXPECT_NEAR_REL_HPP
#define TEST_UNIT_MATH_TEST_EXPECT_NEAR_REL_HPP

#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

namespace stan {
namespace test {

namespace internal {

/**
 * Test that the specified values are within the specified tolerance
 * on relative error, and if not, fail the embedded google test.
 *
 * <p>Relative error is defined to be the error `u - v` rescaled by the
 * average absolute value,
 * `rel_err(u, v) = (u - v) / (0.5 * (abs(u) * + abs(v))).`
 *
 * <p>If at least one of `u` or `v` is zero, the absolute error is
 * tested at the specified tolerance, because the relative error
 * reduces to a constant.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param msg failure message for reporting
 * @param x1 first argument
 * @param x2 second argument
 * @param tol relative tolerance
 */
template <typename T1, typename T2>
void expect_near_rel_finite(const std::string& msg, const T1& x1, const T2& x2,
                            double tol = 1e-8) {
  using stan::math::fabs;
  // if both zero, can just return
  // if only one is zero, test that the non-zero one is close to zero,
  // because general case reduces to 2 if x1 = 0 and x2 != 0 and vice-versa
  if (x1 == 0 && x2 == 0)
    return;
  if (x1 == 0 || x2 == 0) {
    EXPECT_NEAR(x1, x2, tol) << "expect_near_rel_finite(" << x1 << ", " << x2
                             << ", tolerance = " << tol << ")"
                             << "    in: " << msg << std::endl;
    return;
  }
  auto avg = 0.5 * (fabs(x1) + fabs(x2));
  auto relative_diff = (x1 - x2) / avg;
  EXPECT_NEAR(0, relative_diff, tol)
      << "expect_near_rel_finite(" << x1 << ", " << x2
      << ", tolerance = " << tol << ")"
      << "; relative diff = " << relative_diff << std::endl
      << "    in: " << msg << std::endl;
}

template <typename T1, typename T2, int R, int C>
void expect_near_rel_finite(const std::string& msg,
                            const Eigen::Matrix<T1, R, C>& x1,
                            const Eigen::Matrix<T2, R, C>& x2) {
  EXPECT_EQ(x1.rows(), x2.rows());
  EXPECT_EQ(x1.cols(), x2.cols());
  for (int i = 0; i < x1.size(); ++i)
    expect_near_rel_finite(msg, x1(i), x2(i));
}

template <typename T1, typename T2>
void expect_near_rel_finite(const std::string& msg, const std::vector<T1>& x1,
                            const std::vector<T2>& x2) {
  EXPECT_EQ(x1.size(), x2.size());
  for (size_t i = 0; i < x1.size(); ++i)
    expect_near_rel_finite(x1[i], x2[i]);
}

}  // namespace internal

/**
 * Test that scalars x1 and x2 are within the specified relative
 * tolerance, with identity behavior for infinite and NaN values.
 * Relative tolerance is defined in the documentation for
 * `expect_near_rel_finite`.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param msg message indicating what is being tested to print in case
 * of error
 * @param x1 first argument to test
 * @param x2 second argument to test
 * @param tol relative tolerance
 */
template <typename T1, typename T2>
void expect_near_rel(const std::string& msg, const T1& x1, const T2& x2,
                     double tol = 1e-8) {
  if (stan::math::is_nan(x1) || stan::math::is_nan(x2))
    EXPECT_TRUE(stan::math::is_nan(x1) && stan::math::is_nan(x2))
        << "expect_near_rel(" << x1 << ", " << x2 << ")" << std::endl
        << msg << std::endl;
  else if (stan::math::is_inf(x1) || stan::math::is_inf(x2))
    EXPECT_EQ(x1, x2) << "expect_near_rel(" << x1 << ", " << x2 << ")"
                      << std::endl
                      << msg << std::endl;
  else
    internal::expect_near_rel_finite(msg, x1, x2, tol);
}

/**
 * Tests that matrices (or vectors) x1 and x2 are same size and
 * have near values up to specified relative tolerance, accounting for
 * infinite, not-a-number, and zero values as specified in
 * `expect_near_rel_finite`.
 *
 * @tparam T1 type of scalars in first matrix
 * @tparam T2 type of scalars in second matrix
 * @tparam R row specification of matrix
 * @tparam C column specification of the matrix
 * @param msg failure message
 * @param x1 first matrix to test
 * @param x2 second matrix to test
 * @param tol relative tolerance
 */
template <typename T1, typename T2, int R, int C>
void expect_near_rel(const std::string& msg, const Eigen::Matrix<T1, R, C>& x1,
                     const Eigen::Matrix<T2, R, C>& x2, double tol = 1e-8) {
  EXPECT_EQ(x1.rows(), x2.rows()) << "expect_near_rel (Eigen::Matrix)"
                                  << " rows must be same size."
                                  << " x1.rows() = " << x1.rows()
                                  << "; x2.rows() = " << x2.rows() << std::endl
                                  << msg << std::endl;
  EXPECT_EQ(x1.cols(), x2.cols())
      << "expect_near_rel:"
      << "cols must be same size."
      << "x1.cols() = " << x1.cols() << "x2.cols() = " << x2.cols() << ")"
      << std::endl
      << msg << std::endl;
  std::string msg2 = "expect_near_rel; require items x1(i) = x2(i): " + msg;
  for (int i = 0; i < x1.size(); ++i)
    expect_near_rel(msg2, x1(i), x2(i), tol);
}

template <typename T1, typename T2>
void expect_near_rel(const std::string& msg, const std::vector<T1>& x1,
                     const std::vector<T2>& x2, double tol = 1e-8) {
  EXPECT_EQ(x1.size(), x2.size()) << "expect_near_rel (std::vector):"
                                  << " vectors must be same size."
                                  << " x1.size() = " << x1.size()
                                  << "; x2.size() = " << x2.size() << std::endl
                                  << msg << std::endl;
  std::string msg2 = "expect_near_rel; requite items x1[i] = x2[i]: " + msg;
  for (size_t i = 0; i < x1.size(); ++i)
    expect_near_rel(msg2, x1[i], x2[i], tol);
}

}  // namespace test
}  // namespace stan
#endif
