#ifndef STAN_MATH_REV_FUN_RFFT_HPP
#define STAN_MATH_REV_FUN_RFFT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/fft.hpp>
#include <stan/math/prim/fun/rfft.hpp>
#include <stan/math/prim/fun/to_complex.hpp>
#include <stan/math/prim/fun/get_real.hpp>
#include <Eigen/Dense>
#include <complex>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the discrete Fourier transform of the specified real
 * vector.
 *
 * The adjoint computation is given by
 * ```
 * adjoint(x) += length(y) * inv_rfft(adjoint(y))
 * ```
 *
 * If the input is of size zero, the result is a size zero vector.
 *
 * @tparam V type of real vector argument
 * @param[in] x vector to transform
 * @return discrete Fourier transform of `x`
 */
template <typename V, require_eigen_vector_t<V>* = nullptr,
          require_var_t<scalar_type_t<value_type_t<V>>>* = nullptr>
inline Eigen::Matrix<complex_return_t<value_type_t<V>>, -1, 1> rfft(
    const V& x) {
  using return_t = Eigen::Matrix<complex_return_t<value_type_t<V>>, -1, 1>;
  if (unlikely(x.size() <= 1)) {
    return to_complex(x, 0);
  }

  arena_t<V> arena_v = x;
  arena_t<return_t> res = rfft(arena_v.val());

  reverse_pass_callback([arena_v, res]() mutable {
    // NB: NOT inv_rfft since adjoints are not necessarily in the half spectrum
    auto adj_inv_fft = inv_fft(to_complex(res.real().adj(), res.imag().adj()));
    adj_inv_fft *= res.size();
    arena_v.adj() += adj_inv_fft.real();
  });

  return return_t(res);
}

/**
 * Return the real-valued inverse discrete Fourier transform of
 * the specified complex vector.
 *
 * The adjoint computation is given by
 * ```
 * adjoint(y) += (1 / length(x)) * rfft(adjoint(x))
 * ```
 *
 * If the input is of size zero, the result is a size zero vector.
 *
 * @tparam V type of complex vector argument
 * @param[in] y vector to inverse transform
 * @return inverse discrete Fourier transform of `y`
 */
// template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
//           require_var_t<base_type_t<value_type_t<V>>>* = nullptr>
// inline Eigen::Matrix<base_type_t<V>, -1, 1> inv_rfft(const V& y) {
//   using return_t = Eigen::Matrix<base_type_t<V>, -1, 1>;
//   if (unlikely(y.size() <= 1)) {
//     return get_real(y);
//   }

//   arena_t<V> arena_v = y;
//   arena_t<return_t> res
//       = inv_rfft(to_complex(arena_v.real().val(), arena_v.imag().val()));

//   reverse_pass_callback([arena_v, res]() mutable {
//     auto adj_fft = rfft(res.adj());
//     adj_fft /= res.size();

//     arena_v.real().adj() += adj_fft.real();
//     arena_v.imag().adj() += adj_fft.imag();
//   });
//   return return_t(res);
// }

/**
 * Return the two-dimensional discrete Fourier transform of the
 * specified real matrix.  The 2D discrete Fourier transform first
 * runs the discrete Fourier transform on the each row, then on each
 * column of the result.
 *
 * The adjoint computation is given by
 * ```
 * adjoint(x) += size(y) * inv_fft2(adjoint(y))
 * ```
 *
 * @tparam M type of real matrix argument
 * @param[in] x matrix to transform
 * @return discrete 2D Fourier transform of `x`
 */
template <typename M, require_eigen_dense_dynamic_t<M>* = nullptr,
          require_var_t<scalar_type_t<M>>* = nullptr>
inline Eigen::Matrix<complex_return_t<value_type_t<M>>, -1, -1> rfft2(
    const M& x) {
  using return_t = Eigen::Matrix<complex_return_t<value_type_t<M>>, -1, -1>;
  arena_t<M> arena_v = x;
  arena_t<return_t> res = rfft2(arena_v.val());

  reverse_pass_callback([arena_v, res]() mutable {
    // NB: NOT inv_rfft since adjoints are not necessarily in the half spectrum
    auto adj_inv_fft
        = inv_fft2(to_complex(res.real().adj(), res.imag().adj()));
    adj_inv_fft *= res.size();
    arena_v.adj() += adj_inv_fft.real();
  });

  return return_t(res);
}

/**
 * Return the two-dimensional real-valued inverse discrete Fourier
 * transform of the specified complex matrix.  The 2D inverse
 * discrete Fourier transform first runs the 1D inverse Fourier
 * transform on the columns, and then on the resulting rows.
 * The composition of the FFT and inverse FFT (or vice-versa)
 * is the identity.
 *
 * The adjoint computation is given by
 * ```
 * adjoint(y) += (1 / size(x)) * rfft2(adjoint(x))
 * ```
 *
 * @tparam M type of complex matrix argument
 * @param[in] y matrix to inverse trnasform
 * @return inverse discrete 2D Fourier transform of `y`
 */
template <typename M, require_eigen_dense_dynamic_vt<is_complex, M>* = nullptr,
          require_var_t<base_type_t<value_type_t<M>>>* = nullptr>
inline Eigen::Matrix<base_type_t<M>, -1, -1> inv_rfft2(const M& y) {
  using return_t = Eigen::Matrix<base_type_t<M>, -1, -1>;
  arena_t<M> arena_v = y;
  arena_t<return_t> res
      = inv_rfft2(to_complex(arena_v.real().val(), arena_v.imag().val()));

  reverse_pass_callback([arena_v, res]() mutable {
    auto adj_fft = rfft2(res.adj());
    adj_fft /= res.size();

    arena_v.real().adj() += adj_fft.real();
    arena_v.imag().adj() += adj_fft.imag();
  });
  return return_t(res);
}

}  // namespace math
}  // namespace stan
#endif
