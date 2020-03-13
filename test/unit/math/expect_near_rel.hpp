#ifndef TEST_UNIT_MATH_TEST_EXPECT_NEAR_REL_HPP
#define TEST_UNIT_MATH_TEST_EXPECT_NEAR_REL_HPP

#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/relative_tolerance.hpp>
#include <string>
#include <vector>

namespace stan {
namespace test {

namespace internal {

/**
 * Test that the specified values are within the specified tolerance
 * on relative error, and if not, fail the embedded google test.
 *
 * Uses relative_tolerance::inexact to compute the tolerance.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param msg failure message for reporting
 * @param x1 first argument
 * @param x2 second argument
 * @param tol relative tolerance
 */
template <typename T1, typename T2, require_all_stan_scalar_t<T1, T2>...>
void expect_near_rel_finite(const std::string& msg, const T1& x1, const T2& x2,
                            const relative_tolerance tol
                            = relative_tolerance()) {
  double tol_val = tol.inexact(x1, x2);
  EXPECT_NEAR(x1, x2, tol_val)
      << "expect_near_rel_finite in: " << msg << std::endl;
}

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>...>
void expect_near_rel_finite(const std::string& msg, const EigMat1& x1,
                            const EigMat2& x2) {
  EXPECT_EQ(x1.rows(), x2.rows());
  EXPECT_EQ(x1.cols(), x2.cols());
  auto x1_eval = x1.eval();
  auto x2_eval = x2.eval();
  for (int i = 0; i < x1.size(); ++i)
    expect_near_rel_finite(msg, x1_eval(i), x2_eval(i));
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
template <typename T1, typename T2, require_all_stan_scalar_t<T1, T2>...>
void expect_near_rel(const std::string& msg, const T1& x1, const T2& x2,
                     relative_tolerance tol = relative_tolerance()) {
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
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>...>
void expect_near_rel(const std::string& msg, EigMat1&& x1, EigMat2&& x2,
                     relative_tolerance tol = relative_tolerance()) {
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
  auto x1_eval = x1.eval();
  auto x2_eval = x2.eval();
  for (int i = 0; i < x1.size(); ++i)
    expect_near_rel(msg2, x1_eval(i), x2_eval(i), tol);
}

template <typename T1, typename T2>
void expect_near_rel(const std::string& msg, const std::vector<T1>& x1,
                     const std::vector<T2>& x2,
                     relative_tolerance tol = relative_tolerance()) {
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
