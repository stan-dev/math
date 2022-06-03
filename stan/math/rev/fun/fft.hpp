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

template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_var_t<base_type_t<value_type_t<V>>>* = nullptr>
inline plain_type_t<V> fft(const V& v) {
  if (unlikely(v.size() <= 1)) {
    return plain_type_t<V>(v);
  }

  arena_t<V> arena_v = v;
  arena_t<V> res = fft(to_complex(arena_v.real().val(), arena_v.imag().val()));

  reverse_pass_callback([arena_v, res]() mutable {
    // adjoint(x) += length(y) * ifft(adjoint(y))
    auto adj_inv_fft = inv_fft(to_complex(res.real().adj(), res.imag().adj()));
    adj_inv_fft *= res.size();
    arena_v.real().adj() += adj_inv_fft.real();
    arena_v.imag().adj() += adj_inv_fft.imag();
  });

  return plain_type_t<V>(res);
}

template <typename V, require_eigen_vector_vt<is_complex, V>* = nullptr,
          require_var_t<base_type_t<value_type_t<V>>>* = nullptr>
inline plain_type_t<V> inv_fft(const V& v) {
  if (unlikely(v.size() <= 1)) {
    return plain_type_t<V>(v);
  }

  arena_t<V> arena_v = v;
  arena_t<V> res
      = inv_fft(to_complex(arena_v.real().val(), arena_v.imag().val()));

  reverse_pass_callback([arena_v, res]() mutable {
    // adjoint(x) += (1 / length(y)) * fft(adjoint(y))
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
