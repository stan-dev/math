#ifndef STAN_MATH_FWD_FUN_CEIL_HPP
#define STAN_MATH_FWD_FUN_CEIL_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the ceiling of the specified variable.
 *
 * The derivative of the ceiling function is defined and
 * zero everywhere but at integers, and we set them to zero for
 * convenience.
 *
 * @tparam T Inner type of the fvar
 * @param x fvar<T> variable
 * @return Ceiling of the variable.
 */
template <typename T>
inline fvar<T> ceil(const fvar<T>& x) {
  using std::ceil;
  return fvar<T>(ceil(x.val_), 0);
}

/**
 * Return the ceiling of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of fvar
 * @return Ceiling of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto ceil(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        const Eigen::Ref<const T_plain>& v_ref = v;

        plain_type_t<decltype(v)> result(v_ref.rows(), v_ref.cols());
        result.val() = ceil(v_ref.val());
        result.d().setZero();

        return result;
});
}
}  // namespace math
}  // namespace stan
#endif
