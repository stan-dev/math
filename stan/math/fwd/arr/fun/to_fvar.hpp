#ifndef STAN_MATH_FWD_ARR_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_ARR_FUN_TO_FVAR_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/to_fvar.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T>& v) {
  std::vector<fvar<T>> x(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    x[i] = T(v[i]);
  }
  return x;
}

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T>& v,
                                    const std::vector<T>& d) {
  std::vector<fvar<T>> x(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    x[i] = fvar<T>(v[i], d[i]);
  }
  return x;
}

/**
 * Specialization of to_fvar for const fvar input
 *
 * @tparam The inner type of the fvar.
 * @param[in,out] v A vector of forward automatic differentiation variable.
 * @return The input vector of forward automatic differentiation variable.
 */
template <typename T>
inline const std::vector<fvar<T>>& to_fvar(const std::vector<fvar<T>>& v) {
  return v;
}

/**
 * Specialization of to_fvar for non-const fvar input
 *
 * @tparam The inner type of the fvar.
 * @param[in,out] v A vector of forward automatic differentiation variable.
 * @return The input vector of forward automatic differentiation variable.
 */
template <typename T>
inline std::vector<fvar<T>>& to_fvar(std::vector<fvar<T>>& v) {
  return v;
}

}  // namespace math
}  // namespace stan
#endif
