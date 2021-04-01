#ifndef STAN_MATH_PRIM_FUN_RANK_HPP
#define STAN_MATH_PRIM_FUN_RANK_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Return the number of components of v less than v[s].
 *
 * @tparam C container type
 * @param[in] v input vector
 * @param[in] s position in vector
 * @return number of components of v less than v[s].
 * @throw std::out_of_range if s is out of range.
 */
template <typename C, require_container_t<C>* = nullptr>
inline int rank(const C& v, int s) {
  check_range("rank", "v", v.size(), s);
  --s;  // adjust for indexing by one
  return apply_vector_unary<C>::reduce(v, [s](const auto& vec) {
    const auto& vec_ref = to_ref(vec);

    return (vec_ref.array() < vec_ref.coeff(s)).template cast<int>().sum();
  });
}

}  // namespace math
}  // namespace stan

#endif
