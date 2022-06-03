#ifndef STAN_MATH_PRIM_FUN_FFT_HPP
#define STAN_MATH_PRIM_FUN_FFT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <unsupported/Eigen/FFT>
#include <Eigen/Dense>
#include <complex>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the discrete Fourier transform of the specified complex
 * vector.
 *
 * Given an input complex vector `x[0:N-1]` of size `N`, the discrete
 * Fourier transform computes entries of the resulting complex
 * vector `y[0:N-1]` by
 *
 * ```
 * y[n] = SUM_{i < N} x[i] * exp(-n * i * 2 * pi * sqrt(-1) / N)
 * ```
 *
 * If the input is of size zero, the result is a size zero vector.
 *
 * @tparam V type of complex vector argument
 * @param[in] x vector to transform
 * @return discrete Fourier transform of `x`
 */
template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_not_var_t<base_type_t<value_type_t<V>>>* = nullptr>
inline Eigen::Matrix<scalar_type_t<V>, -1, 1> fft(const V& x) {
  // copy because fft() requires Eigen::Matrix type
  Eigen::Matrix<scalar_type_t<V>, -1, 1> xv = x;
  if (xv.size() <= 1)
    return xv;
  Eigen::FFT<base_type_t<V>> fft;
  return fft.fwd(xv);
}

/**
 * Return the inverse discrete Fourier transform of the specified
 * complex vector.
 *
 * Given an input complex vector `y[0:N-1]` of size `N`, the inverse
 * discrete Fourier transform computes entries of the resulting
 * complex vector `x[0:N-1]` by
 *
 * ```
 * x[n] = SUM_{i < N} y[i] * exp(n * i * 2 * pi * sqrt(-1) / N)
 * ```
 *
 * If the input is of size zero, the result is a size zero vector.
 * The only difference between the discrete DFT and its inverse is
 * the sign of the exponent.
 *
 * @tparam V type of complex vector argument
 * @param[in] y vector to inverse transform
 * @return inverse discrete Fourier transform of `y`
 */
template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_not_var_t<base_type_t<value_type_t<V>>>* = nullptr>
inline Eigen::Matrix<scalar_type_t<V>, -1, 1> inv_fft(const V& y) {
  // copy because fft() requires Eigen::Matrix type
  Eigen::Matrix<scalar_type_t<V>, -1, 1> yv = y;
  if (y.size() <= 1)
    return yv;
  Eigen::FFT<base_type_t<V>> fft;
  return fft.inv(yv);
}

/**
 * Return the two-dimensional discrete Fourier transform of the
 * specified complex matrix.  The 2D discrete Fourier transform first
 * runs the discrete Fourier transform on the each row, then on each
 * column of the result.
 *
 * @tparam M type of complex matrix argument
 * @param[in] x matrix to transform
 * @return discrete 2D Fourier transform of `x`
 */
template <typename M, require_eigen_dense_dynamic_vt<is_complex, M>* = nullptr>
inline Eigen::Matrix<scalar_type_t<M>, -1, -1> fft2(const M& x) {
  Eigen::Matrix<scalar_type_t<M>, -1, -1> y(x.rows(), x.cols());
  for (int i = 0; i < y.rows(); ++i)
    y.row(i) = fft(x.row(i));
  for (int j = 0; j < y.cols(); ++j)
    y.col(j) = fft(y.col(j));
  return y;
}

/**
 * Return the two-dimensional inverse discrete Fourier transform of
 * the specified complex matrix.  The 2D inverse discrete Fourier
 * transform first runs the 1D inverse Fourier transform on the
 * columns, and then on the resulting rows.  The composition of the
 * FFT and inverse FFT (or vice-versa) is the identity.
 *
 * @tparam M type of complex matrix argument
 * @param[in] y matrix to inverse trnasform
 * @return inverse discrete 2D Fourier transform of `y`
 */
template <typename M, require_eigen_dense_dynamic_vt<is_complex, M>* = nullptr>
inline Eigen::Matrix<scalar_type_t<M>, -1, -1> inv_fft2(const M& y) {
  Eigen::Matrix<scalar_type_t<M>, -1, -1> x(y.rows(), y.cols());
  for (int j = 0; j < x.cols(); ++j)
    x.col(j) = inv_fft(y.col(j));
  for (int i = 0; i < x.rows(); ++i)
    x.row(i) = inv_fft(x.row(i));
  return x;
}

}  // namespace math
}  // namespace stan

#endif
