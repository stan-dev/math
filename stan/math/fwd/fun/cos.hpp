#ifndef STAN_MATH_FWD_FUN_COS_HPP
#define STAN_MATH_FWD_FUN_COS_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> cos(const fvar<T>& x) {
  using std::cos;
  using std::sin;
  return fvar<T>(cos(x.val_), x.d_ * -sin(x.val_));
}

template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto cos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        const Eigen::Ref<const T_plain>& v_ref = v;
        auto vals = v_ref.val().eval();

        T_plain result(v_ref.rows(), v_ref.cols());
        result.val() = cos(vals);
        result.d().array() = v_ref.d().array() * -sin(vals).array();

        return result;
});
}
}  // namespace math
}  // namespace stan
#endif
