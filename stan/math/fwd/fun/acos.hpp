#ifndef STAN_MATH_FWD_FUN_ACOS_HPP
#define STAN_MATH_FWD_FUN_ACOS_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the arc cosine of a variable, in radians (cmath).
 *
 * The derivative is defined by:
   \f[
   \frac{\partial \, \arccos(x)}{\partial x} = -\frac{1}{\sqrt{1-x^2}}
   \f]
 *
 * @tparam T Inner type of the fvar
 * @param x fvar<T> variable
 * @return Arc cosine of variable, in radians.
 */
template <typename T>
inline fvar<T> acos(const fvar<T>& x) {
  using std::acos;
  using std::sqrt;
  return fvar<T>(acos(x.val_), x.d_ / -sqrt(1 - square(x.val_)));
}

/**
 * Return the arc cosine of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of fvar
 * @return Arc cosine of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto acos(const Container& x) {
  return apply_vector_unary<Container>::apply(x, [](const auto& v) {
    using T_plain = plain_type_t<decltype(v)>;
    const Eigen::Ref<const T_plain>& v_ref = v;
    auto vals = v_ref.val().eval();

    T_plain result(v_ref.rows(), v_ref.cols());
    result.val() = acos(vals);
    result.d().array()
        = v_ref.d().array() * -inv_sqrt(1 - square(vals).array());

    return result;
  });
}
}  // namespace math
}  // namespace stan
#endif
