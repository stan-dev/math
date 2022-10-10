#ifndef STAN_MATH_FWD_FUN_FFT_HPP
#define STAN_MATH_FWD_FUN_FFT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/fft.hpp>
#include <stan/math/prim/fun/to_complex.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <Eigen/Dense>
#include <complex>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<V>>>* = nullptr>
inline plain_type_t<V> fft(const V& x) {
  if (unlikely(x.size() <= 1)) {
    return plain_type_t<V>(x);
  }

  auto res = fft(to_complex(x.real().val(), x.imag().val()));
  auto dres = fft(to_complex(x.real().d(), x.imag().d()));

  return plain_type_t<V>(to_complex(to_fvar(res.real(), dres.real()),
                                    to_fvar(res.imag(), dres.imag())));
}

template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<V>>>* = nullptr>
inline plain_type_t<V> inv_fft(const V& x) {
  if (unlikely(x.size() <= 1)) {
    return plain_type_t<V>(x);
  }

  auto res = inv_fft(to_complex(x.real().val(), x.imag().val()));
  auto dres = inv_fft(to_complex(x.real().d(), x.imag().d()));

  return plain_type_t<V>(to_complex(to_fvar(res.real(), dres.real()),
                                    to_fvar(res.imag(), dres.imag())));
}

template <typename M, require_eigen_dense_dynamic_vt<is_complex, M>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<M>>>* = nullptr>
inline plain_type_t<M> fft2(const M& x) {
  if (unlikely(x.size() <= 1)) {
    return plain_type_t<M>(x);
  }

  auto res = fft2(to_complex(x.real().val(), x.imag().val()));
  auto dres = fft2(to_complex(x.real().d(), x.imag().d()));

  return plain_type_t<M>(to_complex(to_fvar(res.real(), dres.real()),
                                    to_fvar(res.imag(), dres.imag())));
}

template <typename M, require_eigen_dense_dynamic_vt<is_complex, M>* = nullptr,
          require_fvar_t<base_type_t<value_type_t<M>>>* = nullptr>
inline plain_type_t<M> inv_fft2(const M& x) {
  if (unlikely(x.size() <= 1)) {
    return plain_type_t<M>(x);
  }

  auto res = inv_fft2(to_complex(x.real().val(), x.imag().val()));
  auto dres = inv_fft2(to_complex(x.real().d(), x.imag().d()));

  return plain_type_t<M>(to_complex(to_fvar(res.real(), dres.real()),
                                    to_fvar(res.imag(), dres.imag())));
}

}  // namespace math
}  // namespace stan

#endif
