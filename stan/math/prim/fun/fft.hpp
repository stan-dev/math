#ifndef STAN_MATH_PRIM_FUN_FFT_HPP
#define STAN_MATH_PRIM_FUN_FFT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <unsupported/Eigen/FFT>
#include <complex>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

  /**
   * Return the discrete Fourier transform of the specified complex
   * vector.  The input vector may be considered to be in the time
   * domain and the output will be in the frequency domain. 
   *
   * Given an input complex vector `x[0:N-1]` of size `N`, the discrete
   * Fourier transform computes entries of the resulting complex
   * vector `y[0:N-1]` by 
   *
   * ```
   * y[n] = SUM_{i < N} x[i] * exp(n * i * -2 * pi * sqrt(-1) / N)
   * ```
   * 
   * If `y` is of size zero, the result is a size zero vector.
   *
   * @tparam V type of complex vector argument
   * @param x complex time domain vector to transform
   * @return discrete Fourier transform of `x`
   */
  template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr>
  inline Eigen::Matrix<scalar_type_t<V>, -1, 1> fft(const V& x) {
    Eigen::Matrix<scalar_type_t<V>, -1, 1> xv = x;
    if (xv.size() <= 1) return xv;
    Eigen::FFT<typename scalar_type_t<V>::value_type> fft;
    return fft.fwd(xv);
  }


  /**
   * Return the inverse discrete Fourier transform of the specified
   * complex vector.  The input may be considered to be in the
   * frequency domain and the output will be in the time domain.
   *
   * Given an input complex vector `y[0:N-1]` of size `N`, the inverse
   * discrete Fourier transform computes entries of the resulting
   * complex vector `x[0:N-1]` by 
   *
   * ```
   * x[n] = SUM_{i < N} y[i] * exp(n * i * 2 * pi * sqrt(-1) / N)
   * ```
   *
   * The only difference between the discrete DFT and its inverse is
   * the sign of the exponent.
   *
   * If the input is size zero, the output will be size zero.
   *
   * @tparam V type of complex vector argument
   * @param y complex frequency domain vector to inverse transform
   * @return inverse discrete Fourier transform of `x`
   */
  template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr>
  inline Eigen::Matrix<scalar_type_t<V>, -1, 1> inv_fft(const V& y) {
    Eigen::Matrix<scalar_type_t<V>, -1, 1> yv = y;
    if (y.size() <= 1) return yv;
    Eigen::FFT<typename scalar_type_t<V>::value_type> fft;
    return fft.inv(yv);
  }

  /**
   * Return the two-dimensional discrete Fourier transform to the
   * specified complex matrix.  Given a complex matrix (M x N) matrix
   * x, the definition of the 2D discrete Fourier transform is
   *
   * ```
   * y[m, n] = SUM_{i < M, j < N} x[i, j] 
   *                              * exp(i * m * -2 * pi * sqrt(-1) / N)
   *                              * exp(j * n * -2 * pi * sqrt(-1) / M)
   * ```
   *
   * Another way to view the 2D FFT is as first running a 1D FFT on
   * each row, then running a 1D FFT on each resulting column.
   * 
   * @tparam M type of complex matrix argument
   * @param x complex time-domain matrix 
   * @return discrete 2D Fourier transform of `x`
   */
  template <typename M> // , require_eigen_matrix_dynamic_t<is_complex, M>* = nullptr>
  inline Eigen::Matrix<scalar_type_t<M>, -1, -1> fft2(const M& x) {
    if (x.size() <= 1) return x;
    int rows = x.rows();
    int cols = x.cols();
    Eigen::FFT<typename scalar_type_t<M>::value_type> fft;
    Eigen::Matrix<scalar_type_t<M>, -1, -1> y(rows, cols);
    Eigen::Matrix<scalar_type_t<M>, -1, 1> v_temp;
    if (cols > 1) {
      for (int i = 0; i < rows; ++i) {
	// implicit transposition by assignment
	Eigen::Matrix<scalar_type_t<M>, -1, 1> row_i_t = x.row(i); 
	fft.fwd(v_temp, row_i_t);
	y.row(i) = v_temp;
      }
    }
    if (rows > 1) {
      for (int j = 0; j < cols; ++j) {
	Eigen::Matrix<scalar_type_t<M>, -1, 1> col_j = y.col(j);
	fft.fwd(v_temp, col_j);
	y.col(j) = v_temp;
      }
    }
    return y;
  }

  // /**
  //  * Return the two-dimensional discrete Fourier transform to the
  //  * specified complex matrix.  Given a complex matrix (M x N) matrix
  //  * x, the definition of the 2D discrete Fourier transform is
  //  *
  //  * ```
  //  * x[m, n] = SUM_{i < M, j < N} y[i, j] 
  //  *                              * exp(i * m * 2 * pi * sqrt(-1) / N)
  //  *                              * exp(j * n * 2 * pi * sqrt(-1) / M)
  //  * ```
  //  *
  //  * Another way to view the 2D inverse FFT is as first running a 1D
  //  * inverse FFT on each row, then running a 1D invese FFT on each
  //  * resulting column.
  //  * 
  //  * @tparam M type of complex matrix argument
  //  * @param y complex frequency-domain matrix 
  //  * @return inverse discrete 2D Fourier transform of `x`
  //  */
  // template <typename M, require_eigen_matrix_vt<is_complex, M>* = nullptr>
  // inline Eigen::Matrix<scalar_type_t<V>, -1, -1> inv_fft2(const M& y) {
  //   int rows = y.rows();
  //   int cols = y.cols();
  //   Eigen::FFT<typename scalar_type_t<V>::value_type> fft;
  //   Eigen::Matrix<scalar_type_t<V>, -1, -1> x;
  //   Eigen::Matrix<scalar_type_t<V>, -1, 1> v_temp;
  //   for (i = 0; i < rows; ++i) {
  //     fft.inv(v_temp, y.row(i));
  //     x.row(i) = v_temp;
  //   }
  //   for (j = 0; j < cols; ++j) {
  //     fft.inv(v_temp, x.col(j));
  //     x.col(j) = v_temp;
  //   }
  //   return x;
  // }
  
}  // namespace math
}  // namespace stan

#endif
