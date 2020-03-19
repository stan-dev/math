#ifndef STAN_MATH_FWD_FUN_ASIN_HPP
#define STAN_MATH_FWD_FUN_ASIN_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> asin(const fvar<T>& x) {
  using std::asin;
  using std::sqrt;
  return fvar<T>(asin(x.val_), x.d_ / sqrt(1 - square(x.val_)));
}

template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto asin(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        const Eigen::Ref<const T_plain>& v_ref = v;
        auto vals = v_ref.val().eval();

        T_plain result(v_ref.rows(), v_ref.cols());
        result.val() = asin(vals);
        result.d().array() = v_ref.d().array()
                                * inv_sqrt(1-square(vals).array());

        return result;
});
}
}  // namespace math
}  // namespace stan
#endif
