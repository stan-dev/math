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
                            const relative_tolerance tol = relative_tolerance(),
                            const char* x1_name = "x1",
                            const char* x2_name = "x2") {
  double tol_val = tol.inexact(x1, x2);
  EXPECT_NEAR(x1, x2, tol_val)
      << "expect_near_rel_finite in: " << msg << " for " << x1_name << " vs "
      << x2_name << std::endl;
}

template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>...>
void expect_near_rel_finite(const std::string& msg, const EigMat1& x1,
                            const EigMat2& x2, const char* x1_name = "x1",
                            const char* x2_name = "x2") {
  EXPECT_EQ(x1.rows(), x2.rows());
  EXPECT_EQ(x1.cols(), x2.cols());
  auto x1_eval = x1.eval();
  auto x2_eval = x2.eval();
  for (int i = 0; i < x1.size(); ++i) {
    expect_near_rel_finite(msg, x1_eval(i), x2_eval(i), x1_name, x2_name);
  }
}

template <typename T1, typename T2>
void expect_near_rel_finite(const std::string& msg, const std::vector<T1>& x1,
                            const std::vector<T2>& x2,
                            const char* x1_name = "x1",
                            const char* x2_name = "x2") {
  EXPECT_EQ(x1.size(), x2.size());
  for (size_t i = 0; i < x1.size(); ++i) {
    expect_near_rel_finite(x1[i], x2[i], x1_name, x2_name);
  }
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
                     relative_tolerance tol = relative_tolerance(),
                     const char* x1_name = "x1", const char* x2_name = "x2") {
  if (stan::math::is_nan(x1) || stan::math::is_nan(x2)) {
    EXPECT_TRUE(stan::math::is_nan(x1) && stan::math::is_nan(x2))
        << "expect_near_rel(" << x1 << ", " << x2 << ")" << std::endl
        << msg << " for " << x1_name << " vs " << x2_name << std::endl;
  } else if (stan::math::is_inf(x1) || stan::math::is_inf(x2)) {
    EXPECT_EQ(x1, x2) << "expect_near_rel(" << x1 << ", " << x2 << ")"
                      << std::endl
                      << msg << " for " << x1_name << " vs " << x2_name
                      << std::endl;
  } else {
    internal::expect_near_rel_finite(msg, x1, x2, tol, x1_name, x2_name);
  }
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
                     relative_tolerance tol = relative_tolerance(),
                     const char* x1_name = "x1", const char* x2_name = "x2") {
  EXPECT_EQ(x1.rows(), x2.rows())
      << "expect_near_rel (Eigen::Matrix)"
      << " rows must be same size." << x1_name << ".rows() = " << x1.rows()
      << "; " << x2_name << ".rows() = " << x2.rows() << std::endl
      << msg << std::endl;

  EXPECT_EQ(x1.cols(), x2.cols())
      << "expect_near_rel:"
      << "cols must be same size." << x1_name << ".cols() = " << x1.cols()
      << "; " << x2_name << ".cols() = " << x2.cols() << std::endl
      << msg << std::endl;
  auto x1_eval = x1.eval();
  auto x2_eval = x2.eval();
  int sentinal_val = 0;
  for (int j = 0; j < x1.cols(); ++j) {
    for (int i = 0; i < x1.rows(); ++i) {
      std::string msg2 = std::string("expect_near_rel; require items x1(");
      if (stan::is_vector<EigMat1>::value) {
        msg2 += std::to_string(sentinal_val) + ") = x2("
                + std::to_string(sentinal_val) + "): " + msg;
      } else {
        msg2 += std::to_string(i) + ", " + std::to_string(j) + ") = x2("
                + std::to_string(i) + ", " + std::to_string(j) + "): " + msg;
      }
      expect_near_rel(msg2, x1_eval(sentinal_val), x2_eval(sentinal_val), tol,
                      x1_name, x2_name);
      sentinal_val++;
    }
  }
#ifdef STAN_TEST_PRINT_MATRIX_FAILURE
  if (::testing::Test::HasFailure()) {
    Eigen::IOFormat CleanFmt(5, 0, ", ", "\n", "[", "]");
    FAIL() << "\nx1: \n"
           << x1.format(CleanFmt) << "\nx2: \n"
           << x2.format(CleanFmt) << "\n";
  }
#endif
}

/**
 * Tewsts that the elements of the specified standard vectors are
 * relatively near one another by calling `expect_near_rel`
 * recursively.
 *
 * @tparam T1 value type for first vector
 * @tparam T2 value type for second vector
 * @param[in] x1 first vector
 * @param[in] x2 second vector
 * @param[in] tol relative tolerance
 */
template <typename T1, typename T2>
void expect_near_rel(const std::string& msg, const std::vector<T1>& x1,
                     const std::vector<T2>& x2,
                     relative_tolerance tol = relative_tolerance(),
                     const char* x1_name = "x1", const char* x2_name = "x2") {
  EXPECT_EQ(x1.size(), x2.size())
      << "expect_near_rel (std::vector):"
      << " vectors must be same size. " << x1_name << ".size() = " << x1.size()
      << "; " << x2_name << ".size() = " << x2.size() << std::endl
      << msg << std::endl;
  std::string msg2 = std::string("expect_near_rel; require items ") + "x1[i] = "
                     + "x2[i]: " + msg + " for " + x1_name + " vs " + x2_name;
  for (size_t i = 0; i < x1.size(); ++i)
    expect_near_rel(msg2, x1[i], x2[i], tol, x1_name, x2_name);
}

/**
 * Tests that the real and complex parts of the specified complex numbers
 * by calling `expect_near_rel` recursively with the specified message
 * and tolerance.
 *
 * @tparam T1 value type of first complex number
 * @tparam T2 value type of second complex number
 * @param msg[in] message to print under failure
 * @param z1[in] first complex number
 * @param z2[in] second complex number
 * @param tol[in] tolerance for comparison
 */
template <typename T1, typename T2>
void expect_near_rel(const std::string& msg, const std::complex<T1>& z1,
                     const std::complex<T2>& z2,
                     relative_tolerance tol = relative_tolerance(),
                     const char* x1_name = "x1", const char* x2_name = "x2") {
  expect_near_rel(msg, z1.real(), z2.real(), tol, x1_name, x2_name);
  expect_near_rel(msg, z1.imag(), z2.imag(), tol, x1_name, x2_name);
}

/**
 * Tests that the specified real number is near the real part of the
 * specified complex number and the complex number's imaginary part is
 * near zero by calling `expect_near_rel` recursively with the
 * specified message and tolerance.
 *
 * @tparam T1 value type of first number
 * @tparam T2 value type of second complex number
 * @param msg[in] message to print under failure
 * @param x1[in] real number
 * @param z2[in] complex number
 * @param tol[in] tolerance for comparison
 */
template <typename T1, typename T2>
void expect_near_rel(const std::string& msg, const T1& x1,
                     const std::complex<T2>& z2,
                     relative_tolerance tol = relative_tolerance(),
                     const char* x1_name = "x1", const char* x2_name = "x2") {
  expect_near_rel(msg, x1, z2.real(), tol, x1_name, x2_name);
  expect_near_rel(msg, 0, z2.imag(), tol, x1_name, x2_name);
}

/**
 * Tests that the specified real number is near the real part of the
 * specified complex number and the complex number's imaginary part is
 * near zero by calling `expect_near_rel` recursively with the
 * specified message and tolerance.
 *
 * @tparam T1 value type of first number
 * @tparam T2 value type of second complex number
 * @param msg[in] message to print under failure
 * @param z1[in] complex number
 * @param x2[in] real number
 * @param tol[in] tolerance for comparison
 */
template <typename T1, typename T2>
void expect_near_rel(const std::string& msg, const std::complex<T1>& z1,
                     const T2& x2,
                     relative_tolerance tol = relative_tolerance(),
                     const char* x1_name = "x1", const char* x2_name = "x2") {
  expect_near_rel(msg, z1.real(), x2, tol, x1_name, x2_name);
  expect_near_rel(msg, z1.imag(), 0, tol, x1_name, x2_name);
}

}  // namespace test
}  // namespace stan

#define TO_STRING_(x) #x
#define TO_STRING(x) TO_STRING_(x)
#define EXPECT_NEAR_REL(a, b)  \
  stan::test::expect_near_rel( \
      "Error in file: " __FILE__ ", on line: " TO_STRING(__LINE__), a, b);

#endif
