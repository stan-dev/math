#ifndef STAN_MATH_PRIM_FUN_NORM2_HPP
#define STAN_MATH_PRIM_FUN_NORM2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns L2 norm of a vector. For vectors that equals the square-root of the
 * sum of squares of the elements.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase)
 * @param v Vector.
 * @return L2 norm of v.
 */
template <typename Container, require_st_arithmetic<Container>* = nullptr,
          require_container_t<Container>* = nullptr>
inline auto norm2(const Container& x) {
  return apply_vector_unary<ref_type_t<Container>>::reduce(
      to_ref(x), [](const auto& v) { return v.template lpNorm<2>(); });
}

}  // namespace math
}  // namespace stan

#endif
