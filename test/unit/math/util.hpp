#ifndef TEST_UNIT_MATH_UTIL_HPP
#define TEST_UNIT_MATH_UTIL_HPP

#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <complex>
#include <string>
#include <vector>

using mvar_t = Eigen::Matrix<stan::math::var, -1, -1>;

using var_t = stan::math::var;
using fvar_d_t = stan::math::fvar<double>;
using fvar_fvar_d_t = stan::math::fvar<fvar_d_t>;
using fvar_v_t = stan::math::fvar<stan::math::var>;
using fvar_fvar_v_t = stan::math::fvar<fvar_v_t>;

using cdouble_t = std::complex<double>;
using cvar_t = std::complex<stan::math::var>;
using cfvar_d_t = std::complex<fvar_d_t>;
using cfvar_fvar_d_t = std::complex<fvar_fvar_d_t>;
using cfvar_v_t = std::complex<fvar_v_t>;
using cfvar_fvar_v_t = std::complex<fvar_fvar_v_t>;

namespace stan {
namespace test {

/**
 * Return the Eigen vector with the same size and elements as the
 * specified standard vector.  Elements are copied from the specifed
 * input vector.
 *
 * @tparam T type of scalars in containers
 * @param x standard vector to copy
 * @return Eigen vector corresponding to specified standard vector
 */
template <typename T>
Eigen::Matrix<T, -1, 1> to_eigen_vector(const std::vector<T>& x) {
  return Eigen::Map<const Eigen::Matrix<T, -1, 1>>(x.data(), x.size());
}

/**
 * Return the standard vector with the same size and elements as the
 * specified Eigen matrix, vector, or row vector.
 *
 * @tparam T type of scalars in containers
 * @tparam R row specification of input matrix
 * @tparam C column specification of input matrix
 * @param x Eigen matrix, vector, or row vector to copy
 * @return standard vector corresponding to input
 */
template <typename T, int R, int C>
std::vector<T> to_std_vector(const Eigen::Matrix<T, R, C>& x) {
  std::vector<T> y;
  y.reserve(x.size());
  for (int i = 0; i < x.size(); ++i)
    y.push_back(x(i));
  return y;
}

/**
 * Succeed if the specified value is the identity matrix to within the
 * specified tolerance.
 *
 * @tparam T value type of the matrix
 * @pram[in] I matrix to test
 * @param[in] tol absolute tolerance
 */
template <typename T>
void expect_identity_matrix(const T& I, double tol = 1e-8) {
  using stan::math::value_of_rec;
  for (int i = 0; i < I.rows(); ++i) {
    EXPECT_NEAR(1, value_of_rec(I(i, i)), tol);
    for (int j = 0; j < i; ++j) {
      EXPECT_NEAR(0, value_of_rec(I(i, j)), tol);
      EXPECT_NEAR(0, value_of_rec(I(j, i)), tol);
    }
  }
}

/**
 * Succed if the specified real and imaginary values match those of
 * the specified complex value.
 *
 * @param[in] re expected real value
 * @param[in] im expected imaginary value
 * @param[in] z imaginary number to test
 */
template <typename T>
void expect_complex(double re, double im, const std::complex<T>& z) {
  using stan::math::value_of_rec;
  stan::test::expect_near_rel("complex real value", re, value_of_rec(z.real()));
  stan::test::expect_near_rel("complex imag value", im, value_of_rec(z.imag()));
}

/**
 * Succeed if the specified complex value has real and complex types
 * matching those of the expected complex value.
 *
 * @param[in] x expected complex value
 * @param[in] y complex value to test
 */
template <typename T>
void expect_complex(const cdouble_t& x, const std::complex<T>& y) {
  expect_complex(x.real(), x.imag(), y);
}

}  // namespace test
}  // namespace stan
#endif
