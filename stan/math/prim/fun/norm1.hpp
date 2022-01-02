#ifndef STAN_MATH_PRIM_FUN_NORM1_HPP
#define STAN_MATH_PRIM_FUN_NORM1_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns L1 norm of a vector. For vectors that equals the
 * sum of magnitudes of its individual elements.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase)
 * @param v Vector.
 * @return L1 norm of v.
 */
template <typename Container, require_st_arithmetic<Container>* = nullptr,
          require_container_t<Container>* = nullptr>
inline auto norm1(const Container& x) {
  return apply_vector_unary<ref_type_t<Container>>::reduce(
      to_ref(x), [](const auto& v) { return v.template lpNorm<1>(); });
}

}  // namespace math
}  // namespace stan

#endif
