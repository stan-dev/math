#ifndef STAN_MATH_REV_FUN_FFT_HPP
#define STAN_MATH_REV_FUN_FFT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/fft.hpp>
#include <stan/math/prim/fun/to_complex.hpp>
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
 * The adjoint computation is given by
 * ```
 * adjoint(x) += length(y) * inv_fft(adjoint(y))
 * ```
 *
 * If the input is of size zero, the result is a size zero vector.
 *
 * @tparam V type of complex vector argument
 * @param[in] x vector to transform
 * @return discrete Fourier transform of `x`
 */
template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_var_t<base_type_t<value_type_t<V>>>* = nullptr>
inline plain_type_t<V> fft(const V& x) {
  if (unlikely(x.size() <= 1)) {
    return plain_type_t<V>(x);
  }

  arena_t<V> arena_v = x;
  arena_t<V> res = fft(to_complex(arena_v.real().val(), arena_v.imag().val()));

  reverse_pass_callback([arena_v, res]() mutable {
    auto adj_inv_fft = inv_fft(to_complex(res.real().adj(), res.imag().adj()));
    adj_inv_fft *= res.size();
    arena_v.real().adj() += adj_inv_fft.real();
    arena_v.imag().adj() += adj_inv_fft.imag();
  });

  return plain_type_t<V>(res);
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
 *  * The adjoint computation is given by
 * ```
 * adjoint(y) += (1 / length(x)) * fft(adjoint(x))
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
          require_var_t<base_type_t<value_type_t<V>>>* = nullptr>
inline plain_type_t<V> inv_fft(const V& y) {
  if (unlikely(y.size() <= 1)) {
    return plain_type_t<V>(y);
  }

  arena_t<V> arena_v = y;
  arena_t<V> res
      = inv_fft(to_complex(arena_v.real().val(), arena_v.imag().val()));

  reverse_pass_callback([arena_v, res]() mutable {
    auto adj_fft = fft(to_complex(res.real().adj(), res.imag().adj()));
    adj_fft /= res.size();

    arena_v.real().adj() += adj_fft.real();
    arena_v.imag().adj() += adj_fft.imag();
  });
  return plain_type_t<V>(res);
}

}  // namespace math
}  // namespace stan
#endif
