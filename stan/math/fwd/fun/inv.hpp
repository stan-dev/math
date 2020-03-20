#ifndef STAN_MATH_FWD_FUN_INV_HPP
#define STAN_MATH_FWD_FUN_INV_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {
/**
 * Return the inverse of the specified variable (cmath).
 *
 * The derivative is defined by
 *
   \f[
   \frac{\partial\, \mbox{inv}(x)}{\partial x} =
   \begin{cases}
     -\frac{1}{x^2} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T Inner type of the fvar
 * @param x fvar<T> variable
 * @return inverse of variable.
 */
template <typename T>
inline fvar<T> inv(const fvar<T>& x) {
  return fvar<T>(1 / x.val_, -x.d_ / square(x.val_));
}

/**
 * Return the inverse of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of fvar
 * @return Inverse of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto inv(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [](const auto& v) {
    using T_plain = plain_type_t<decltype(v)>;
    const Eigen::Ref<const T_plain>& v_ref = v;
    auto vals = v_ref.val().eval();

    T_plain result(v_ref.rows(), v_ref.cols());
    result.val() = inv(vals);
    result.d().array() = v_ref.d().array() * -square(result.val()).array();

    return result;
  });
}
}  // namespace math
}  // namespace stan
#endif
