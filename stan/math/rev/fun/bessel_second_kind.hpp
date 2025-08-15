#ifndef STAN_MATH_REV_FUN_BESSEL_SECOND_KIND_HPP
#define STAN_MATH_REV_FUN_BESSEL_SECOND_KIND_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/bessel_second_kind.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, require_integral_t<T1>* = nullptr,
          require_var_t<T2>* = nullptr, require_stan_scalar_t<T2>* = nullptr>
inline var bessel_second_kind(T1&& v, T2&& a) {
  double ret_val = bessel_second_kind(v, a.val());
  auto precomp_bessel
      = v * ret_val / a.val() - bessel_second_kind(v + 1, a.val());
  return make_callback_var(ret_val, [precomp_bessel, a](auto& vi) mutable {
    a.adj() += vi.adj() * precomp_bessel;
  });
}

/**
 * Overload with `var_value<Matrix>` for `int`, `std::vector<int>`, and
 * `std::vector<std::vector<int>>`
 */
template <typename T1, typename T2, require_st_integral<T1>* = nullptr,
          require_var_matrix_t<T2>* = nullptr>
inline auto bessel_second_kind(T1&& v, T2&& a) {
  auto ret_val = bessel_second_kind(v, a.val()).array().eval();
  auto v_map = as_array_or_scalar(v);
  auto precomp_bessel
      = to_arena(v_map * ret_val / a.val().array()
                 - bessel_second_kind(v_map + 1, a.val().array()));
  return make_callback_var(
      ret_val.matrix(), [precomp_bessel, a](const auto& vi) mutable {
        a.adj().array() += vi.adj().array() * precomp_bessel;
      });
}

}  // namespace math
}  // namespace stan
#endif
