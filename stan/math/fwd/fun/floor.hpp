#ifndef STAN_MATH_FWD_FUN_FLOOR_HPP
#define STAN_MATH_FWD_FUN_FLOOR_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> floor(const fvar<T>& x) {
  using std::floor;
  return fvar<T>(floor(x.val_), 0);
}

template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto floor(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        const Eigen::Ref<const T_plain>& v_ref = v;

        plain_type_t<decltype(v)> result(v_ref.rows(), v_ref.cols());
        result.val() = floor(v_ref.val());
        result.d().setZero();

        return result;
});
}
}  // namespace math
}  // namespace stan
#endif
