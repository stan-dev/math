#ifndef STAN_MATH_FWD_FUN_FABS_HPP
#define STAN_MATH_FWD_FUN_FABS_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {
/**
 * Return the absolute value of the variable.
 *
 * The derivative is defined by
 *
   \f[
   \frac{\partial\, \mbox{fabs}(x)}{\partial x} =
   \begin{cases}
     -1 & \mbox{if } x < 0 \\
     0 & \mbox{if } x = 0 \\
     1 & \mbox{if } x > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T Inner type of the fvar
 * @param x fvar<T> variable
 * @return Absolute value of variable.
 */
template <typename T>
inline fvar<T> fabs(const fvar<T>& x) {
  using std::fabs;

  if (unlikely(is_nan(value_of(x.val_)))) {
    return fvar<T>(fabs(x.val_), NOT_A_NUMBER);
  } else if (x.val_ > 0.0) {
    return x;
  } else if (x.val_ < 0.0) {
    return fvar<T>(-x.val_, -x.d_);
  } else {
    return fvar<T>(0, 0);
  }
}

/**
 * Return the absolute value of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of fvar
 * @return Absolute value of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto fabs(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [](const auto& v) {
    using T_plain = plain_type_t<decltype(v)>;
    const Eigen::Ref<const T_plain>& v_ref = v;
    auto vals = v_ref.val().eval();

    T_plain result(v_ref.rows(), v_ref.cols());
    result.val() = fabs(vals);
    result.d().array() = vals.array().isNaN().select(
        NOT_A_NUMBER,
        (vals.array() < 0).select(-v_ref.d().array(), v_ref.d().array()));
    return result;
  });
}
}  // namespace math
}  // namespace stan
#endif
