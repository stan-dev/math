#ifndef STAN_MATH_PRIM_FUN_LMULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_LMULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the first argument times the log of the second argument. The
 * result is 0 if both arguments are 0.  The funcgtion is defined by
 * `lmultiply(x, y) = x * log(y)` if `x` or `y` is non-zero and
 * `lmultiply(0, 0) = 0` otherwise.
 *
 * @tparam Scalar1 type of the first argument
 * @tparam Scalar2 type of the second argument
 * @param a first argument
 * @param b second argument
 * @return the first argument times the log of the second argument
 */
template <typename Scalar1, typename Scalar2, require_all_arithmetic_t<Scalar1, Scalar2>* = nullptr>
inline return_type_t<Scalar1, Scalar2> lmultiply(const Scalar1 a, const Scalar2 b) {
  using std::log;
  if (a == 0 && b == 0) {
    return 0;
  }
  return a * log(b);
}

/**
 * Return the result of applying `lmultiply` to the arguments
 * elementwise, with broadcasting if one of the arguments is a scalar.
 * At least one of the arguments must be a container.
 *
 * @tparam Container1 type of the first argument
 * @tparam Container2 type of the second argument
 * @param a first argument
 * @param b second argument
 * @return result of applying `lmultiply` to the arguments
 */
template <typename Container1, typename Container2, require_any_container_t<Container1, Container2>* = nullptr,
          require_all_not_var_matrix_t<Container1, Container2>* = nullptr>
inline auto lmultiply(const Container1& a, const Container2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return lmultiply(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
