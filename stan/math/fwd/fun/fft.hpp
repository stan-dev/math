#ifndef STAN_MATH_FWD_FUN_FFT_HPP
#define STAN_MATH_FWD_FUN_FFT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/fft.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the discrete Fourier transform of the specified complex
 * vector for forward-mode autodiff.
 *
 * @tparam V type of complex vector argument
 * @param[in] x vector to transform
 * @return discrete Fourier transform of `x`
 */
template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<V>>>* = nullptr>
inline Eigen::Matrix<scalar_type_t<V>, -1, 1> fft(V&& x) {
  using scalar_t = scalar_type_t<V>;
  using fvar_t = base_type_t<scalar_t>;
  using complex_t = std::complex<partials_type_t<fvar_t>>;
  decltype(auto) x_ref = to_ref(std::forward<V>(x));
  if (x_ref.size() <= 1) {
    return Eigen::Matrix<scalar_type_t<V>, -1, 1>(x_ref);
  }

  Eigen::Matrix<complex_t, -1, 1> x_val = value_of(x_ref);
  Eigen::Matrix<complex_t, -1, 1> x_d = x_ref.unaryExpr(
      [](const auto& z) { return complex_t(z.real().d(), z.imag().d()); });

  auto y_val = fft(std::move(x_val));
  auto y_d = fft(std::move(x_d));

  using out_t = Eigen::Matrix<scalar_type_t<V>, -1, 1>;
  out_t y
      = y_val.binaryExpr(y_d, [](const complex_t& val, const complex_t& der) {
          return std::complex<fvar_t>{fvar_t(val.real(), der.real()),
                                      fvar_t(val.imag(), der.imag())};
        });
  return y;
}

/**
 * Return the inverse discrete Fourier transform of the specified
 * complex vector for forward-mode autodiff.
 *
 * @tparam V type of complex vector argument
 * @param[in] y vector to inverse transform
 * @return inverse discrete Fourier transform of `y`
 */
template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<V>>>* = nullptr>
inline Eigen::Matrix<scalar_type_t<V>, -1, 1> inv_fft(V&& y) {
  using scalar_t = scalar_type_t<V>;
  using fvar_t = base_type_t<scalar_t>;
  using complex_t = std::complex<partials_type_t<fvar_t>>;
  decltype(auto) y_ref = to_ref(std::forward<V>(y));
  if (y_ref.size() <= 1) {
    return Eigen::Matrix<scalar_type_t<V>, -1, 1>(y_ref);
  }

  Eigen::Matrix<complex_t, -1, 1> y_val = value_of(y_ref);
  Eigen::Matrix<complex_t, -1, 1> y_d = y_ref.unaryExpr(
      [](const auto& z) { return complex_t(z.real().d(), z.imag().d()); });

  auto x_val = inv_fft(std::move(y_val));
  auto x_d = inv_fft(std::move(y_d));

  using out_t = Eigen::Matrix<scalar_type_t<V>, -1, 1>;
  out_t x
      = x_val.binaryExpr(x_d, [](const complex_t& val, const complex_t& der) {
          return std::complex<fvar_t>{fvar_t(val.real(), der.real()),
                                      fvar_t(val.imag(), der.imag())};
        });
  return x;
}

/**
 * Return the two-dimensional discrete Fourier transform of the
 * specified complex matrix for forward-mode autodiff.
 *
 * @tparam M type of complex matrix argument
 * @param[in] x matrix to transform
 * @return discrete 2D Fourier transform of `x`
 */
template <typename M, require_eigen_dense_dynamic_vt<is_complex, M>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<M>>>* = nullptr>
inline Eigen::Matrix<scalar_type_t<M>, -1, -1> fft2(M&& x) {
  using scalar_t = scalar_type_t<M>;
  using fvar_t = base_type_t<scalar_t>;
  using complex_t = std::complex<partials_type_t<fvar_t>>;
  decltype(auto) x_ref = to_ref(std::forward<M>(x));
  Eigen::Matrix<complex_t, -1, -1> x_val = value_of(x_ref);
  Eigen::Matrix<complex_t, -1, -1> x_d = x_ref.unaryExpr(
      [](const auto& z) { return complex_t(z.real().d(), z.imag().d()); });

  auto y_val = fft2(std::move(x_val));
  auto y_d = fft2(std::move(x_d));

  using out_t = Eigen::Matrix<scalar_type_t<M>, -1, -1>;
  out_t y
      = y_val.binaryExpr(y_d, [](const complex_t& val, const complex_t& der) {
          return std::complex<fvar_t>{fvar_t(val.real(), der.real()),
                                      fvar_t(val.imag(), der.imag())};
        });
  return y;
}

/**
 * Return the two-dimensional inverse discrete Fourier transform of
 * the specified complex matrix for forward-mode autodiff.
 *
 * @tparam M type of complex matrix argument
 * @param[in] y matrix to inverse transform
 * @return inverse discrete 2D Fourier transform of `y`
 */
template <typename M, require_eigen_dense_dynamic_vt<is_complex, M>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<M>>>* = nullptr>
inline Eigen::Matrix<scalar_type_t<M>, -1, -1> inv_fft2(M&& y) {
  using scalar_t = scalar_type_t<M>;
  using fvar_t = base_type_t<scalar_t>;
  using complex_t = std::complex<partials_type_t<fvar_t>>;
  decltype(auto) y_ref = to_ref(std::forward<M>(y));
  Eigen::Matrix<complex_t, -1, -1> y_val = value_of(y_ref);
  Eigen::Matrix<complex_t, -1, -1> y_d = y_ref.unaryExpr(
      [](const auto& z) { return complex_t(z.real().d(), z.imag().d()); });

  auto x_val = inv_fft2(std::move(y_val));
  auto x_d = inv_fft2(std::move(y_d));

  using out_t = Eigen::Matrix<scalar_type_t<M>, -1, -1>;
  out_t x
      = x_val.binaryExpr(x_d, [](const complex_t& val, const complex_t& der) {
          return std::complex<fvar_t>{fvar_t(val.real(), der.real()),
                                      fvar_t(val.imag(), der.imag())};
        });
  return x;
}

}  // namespace math
}  // namespace stan

#endif
