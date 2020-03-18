#ifndef STAN_MATH_FWD_FUN_EXP_HPP
#define STAN_MATH_FWD_FUN_EXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
namespace math {
template <typename T>
inline fvar<T> exp(const fvar<T>& x) {
  using std::exp;
  return fvar<T>(exp(x.val_), x.d_ * exp(x.val_));
}

template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto exp(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_vec = decltype(v);
        auto exp_val = exp(v.val()).eval();

        return to_fvar(exp_val,exp_val.cwiseProduct(v.d()));
});
}

}  // namespace math
}  // namespace stan
#endif
