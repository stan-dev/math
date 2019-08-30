#ifndef STAN_MATH_FWD_ARR_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_ARR_FUN_TO_FVAR_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/to_fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <utility>
#include <vector>

namespace stan {
namespace math {

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T>& x) {
  std::vector<fvar<T>> x_ret(x.size());
  for (size_t i = 0; i < x.size(); ++i)
    x_ret[i] = T(x[i]);
  return x_ret;
}

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T>& v,
                                    const std::vector<T>& d) {
  std::vector<fvar<T>> x(v.size());
  for (size_t i = 0; i < v.size(); ++i)
    x[i] = fvar<T>(v[i], d[i]);
  return x;
}

/**
 * Specialization of to_fvar for fvar input
 *
 * @tparam The inner type of the fvar.
 * @param[in,out] x A vector of forward automatic differentiation variable.
 * @return The input vector of forward automatic differentiation variable.
 */
template <typename T, enable_if_vector<T>...,
          enable_if_fvar<scalar_type_decay_t<T>>...>
inline auto&& to_fvar(T&& x) {
  return std::forward<T>(x);
}

}  // namespace math
}  // namespace stan
#endif
