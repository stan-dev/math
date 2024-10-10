#ifndef STAN_MATH_REV_CORE_ADJOINT_HPP
#define STAN_MATH_REV_CORE_ADJOINT_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <stan/math/rev/meta/is_var_complex.hpp>
#include <stan/math/rev/meta/is_vari.hpp>

namespace stan {
namespace math {
template <typename Acc, typename T>
inline constexpr auto adjoint(T&& x) noexcept {
  if constexpr (is_var_complex_v<Acc>) {
    return x.adj();
  } else if constexpr (!is_var_complex_v<Acc>) {
    if constexpr (is_complex_v<typename std::decay_t<T>::value_type>) {
      return x.adj().real();
    } else if constexpr (std::is_floating_point_v<typename std::decay_t<T>::value_type>) {
      return x.adj();
    } else {
      static_assert(1, "NOOOO");
    }
  }
}
}
}

#endif
