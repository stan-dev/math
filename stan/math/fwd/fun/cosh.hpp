#ifndef STAN_MATH_FWD_FUN_COSH_HPP
#define STAN_MATH_FWD_FUN_COSH_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> cosh(const fvar<T>& x) {
  using std::cosh;
  using std::sinh;
  return fvar<T>(cosh(x.val_), x.d_ * sinh(x.val_));
}

template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto cosh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        const Eigen::Ref<const T_plain>& v_ref = v;
        auto vals = v_ref.val().eval();

        T_plain result(v_ref.rows(), v_ref.cols());
        result.val() = cosh(vals);
        result.d().array() = v_ref.d().array() * sinh(vals).array();

        return result;
});
}

}  // namespace math
}  // namespace stan
#endif
