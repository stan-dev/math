#ifndef STAN_MATH_PRIM_FUN_ANY_HPP
#define STAN_MATH_PRIM_FUN_ANY_HPP

#include <stan/math/prim/meta.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Return true if any values in the input are true.
 *
 * Overload for a single boolean input
 *
 * @param x boolean input
 * @return The input unchanged
 */
inline bool any(bool x) { return x; }

/**
 * Return true if any values in the input are true.
 *
 * Overload for a std::vector<bool> input. The Eigen::Map/apply_vector_unary
 * approach cannot be used as std::vector<bool> types do not have a .data()
 * member and are not always stored contiguously.
 *
 * @param x std::vector of boolean inputs
 * @return Boolean indicating whether any elements are true
 */
inline bool any(const std::vector<bool>& x) {
  return std::any_of(x.begin(), x.end(), [](bool i) { return i; });
}

/**
 * Return true if any values in the input are true.
 *
 * Overload for Eigen types
 *
 * @tparam Eigen type of the input
 * @param x Eigen object of boolean inputs
 * @return Boolean indicating whether any elements are true
 */
template <typename ContainerT,
          require_eigen_st<std::is_integral, ContainerT>* = nullptr>
inline bool any(const ContainerT& x) {
  return x.any();
}

/**
 * Return true if any values in the input are true.
 *
 * Overload for nested containers
 *
 * @tparam Type of container within std::vector
 * @param x Nested container of boolean inputs
 * @return Boolean indicating whether any elements are true
 */
template <typename ContainerT, require_container_t<ContainerT>* = nullptr>
inline bool any(const std::vector<ContainerT>& x) {
  return std::any_of(x.begin(), x.end(), [](const auto& i) { return any(i); });
}

}  // namespace math
}  // namespace stan

#endif
