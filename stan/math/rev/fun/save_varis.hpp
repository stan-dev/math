#ifndef STAN_MATH_REV_FUN_SAVE_VARIS_HPP
#define STAN_MATH_REV_FUN_SAVE_VARIS_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

vari** save_varis(vari** dest) { return dest; }

template <typename... Pargs>
vari** save_varis(vari** dest, const var& x, const Pargs&... args);

template <typename... Pargs>
vari** save_varis(vari** dest, const std::vector<var>& x, const Pargs&... args);

template <typename T, require_t<is_var<scalar_type_t<T>>>..., typename... Pargs>
vari** save_varis(vari** dest, const std::vector<T>& x, const Pargs&... args);

template <typename... Pargs, typename Mat, require_eigen_vt<is_var, Mat>...>
vari** save_varis(vari** dest, const Mat& x, const Pargs&... args);

template <typename R, require_arithmetic_t<scalar_type_t<R>>...,
          typename... Pargs>
vari** save_varis(vari** dest, const R& x, const Pargs&... args);

template <typename... Pargs>
vari** save_varis(vari** dest, const var& x, const Pargs&... args) {
  *dest = x.vi_;
  return save_varis(dest + 1, args...);
}

template <typename... Pargs>
vari** save_varis(vari** dest, const std::vector<var>& x,
                  const Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest[i] = x[i].vi_;
  }
  return save_varis(dest + x.size(), args...);
}

template <typename T, require_t<is_var<scalar_type_t<T>>>..., typename... Pargs>
vari** save_varis(vari** dest, const std::vector<T>& x, const Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest = save_varis(dest, x[i]);
  }
  return save_varis(dest, args...);
}

template <typename... Pargs, typename Mat, require_eigen_vt<is_var, Mat>...>
vari** save_varis(vari** dest, const Mat& x, const Pargs&... args) {
  for (size_t i = 0; i < x.size(); ++i) {
    dest[i] = x(i).vi_;
  }
  return save_varis(dest + x.size(), args...);
}

template <typename R, require_arithmetic_t<scalar_type_t<R>>...,
          typename... Pargs>
vari** save_varis(vari** dest, const R& x, const Pargs&... args) {
  return save_varis(dest, args...);
}

}  // namespace internal

}  // namespace math
}  // namespace stan
#endif
