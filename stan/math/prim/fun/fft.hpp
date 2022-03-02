#ifndef STAN_MATH_PRIM_FUN_FFT_HPP
#define STAN_MATH_PRIM_FUN_FFT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <unsupported/Eigen/FFT>
#include <complex>
#include <vector>

namespace stan {
namespace math {

  /**
   * Return the forward discrete Fourier transform of the specified
   * complex vector.  If the input is size zero, the output will be
   * size zero.
   *
   * Given an input complex vector `x[0:N-1]` of size `N`, the forward
   * FFT computes entries of the resulting complex vector `y[0:N-1]` by
   *
   * ```
   * y[i] = SUM_{j < N} x[j] exp(-2 * pi * i * j * sqrt(-1) / n)
   * ```
   *
   * @tparam T scalar type of real and imaginary components
   * @param x complex time vector to transform
   * @return discrete Fourier transform of `x`
   */
  template <typename T>
  Eigen::Matrix<std::complex<T>, -1, 1>
  fft(const Eigen::Matrix<std::complex<T>, -1, 1>& x) {
    if (x.size() <= 1) return x;
    Eigen::FFT<T> fft;
    auto y = fft.fwd(x);
    return y;
  }


  /**
   * Return the inverse discrete Fourier transform of the specified
   * complex vector.  If the input is size zero, the output will be
   * size zero.
   *
   * Given an input complex vector `y[0:N-1]` of size `N`, the reverse
   * FFT computes entries of the resulting complex vector `x[0:N-1]` by
   *
   * ```
   * x[i] = SUM_{j < N} y[j] exp(2 * pi * i * j * sqrt(-1) / n)
   * ```
   *
   * @tparam T scalar type of real and imaginary components
   * @param y complex frequency domain vector to inverse transform
   * @return inverse discrete Fourier transform of `x`
   */
  template <typename T>
  Eigen::Matrix<std::complex<T>, -1, 1>
  inv_fft(const Eigen::Matrix<std::complex<T>, -1, 1>& y) {
    if (y.size() <= 1) return y;
    Eigen::FFT<T> fft;
    auto x = fft.inv(y);
    return x;
  }
  
}
}

#endif
