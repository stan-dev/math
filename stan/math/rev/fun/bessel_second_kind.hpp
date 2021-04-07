#ifndef STAN_MATH_REV_FUN_BESSEL_SECOND_KIND_HPP
#define STAN_MATH_REV_FUN_BESSEL_SECOND_KIND_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/bessel_second_kind.hpp>

namespace stan {
namespace math {

inline var bessel_second_kind(int v, const var& a) {
  double ret_val = bessel_second_kind(v, a.val());
  auto precomp_bessel
      = v * ret_val / a.val() - bessel_second_kind(v + 1, a.val());
  return make_callback_var(ret_val, [precomp_bessel, a](auto& vi) mutable {
    a.adj() += vi.adj() * precomp_bessel;
  });
}

template <typename T, require_eigen_t<T>* = nullptr>
inline auto bessel_second_kind(int v, const var_value<T>& a) {
  auto ret_val = bessel_second_kind(v, a.val()).array().eval();
  auto precomp_bessel = to_arena(v * ret_val / a.val().array()
                                 - bessel_second_kind(v + 1, a.val()).array());
  return make_callback_var(
      ret_val.matrix(), [precomp_bessel, a](const auto& vi) mutable {
        a.adj().array() += vi.adj().array() * precomp_bessel;
      });
}

template <typename T, require_eigen_t<T>* = nullptr>
inline auto bessel_second_kind(const std::vector<int>& v,
                               const var_value<T>& a) {
  auto ret_val = bessel_second_kind(v, a.val()).array().eval();
  Eigen::Map<const Eigen::Array<int, -1, 1>> v_map(v.data(), v.size());
  auto precomp_bessel
      = to_arena(v_map.template cast<double>() * ret_val / a.val().array()
                 - bessel_second_kind(v_map + 1, a.val().array()));
  return make_callback_var(
      ret_val.matrix(), [precomp_bessel, a](const auto& vi) mutable {
        a.adj().array() += vi.adj().array() * precomp_bessel;
      });
}

}  // namespace math
}  // namespace stan
#endif
