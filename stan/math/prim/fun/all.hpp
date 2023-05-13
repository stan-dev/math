#ifndef STAN_MATH_PRIM_FUN_ALL_HPP
#define STAN_MATH_PRIM_FUN_ALL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Return true if all values in the input are true.
 *
 * Overload for a single boolean input
 *
 * @param x boolean input
 * @return The input unchanged
 */
inline bool all(bool x) { return x; }

/**
 * Return true if all values in the input are true.
 *
 * Overload for Eigen types
 *
 * @tparam Eigen type of the input
 * @param x Eigen object of boolean inputs
 * @return Boolean indicating whether all elements are true
 */
template <typename ContainerT,
          require_eigen_st<std::is_integral, ContainerT>* = nullptr>
inline bool all(const ContainerT& x) {
  return x.all();
}

/**
 * Return true if all values in the input are true.
 *
 * Overload for a std::vector/nested inputs. The Eigen::Map/apply_vector_unary
 * approach cannot be used as std::vector<bool> types do not have a .data()
 * member and are not always stored contiguously.
 *
 * @tparam Type of container
 * @param x Nested container of boolean inputs
 * @return Boolean indicating whether all elements are true
 */
template <typename ContainerT,
          require_std_vector_st<std::is_integral, ContainerT>* = nullptr>
inline bool all(const ContainerT& x) {
  return std::all_of(x.begin(), x.end(), [](const auto& i) { return all(i); });
}

/**
 * Return true if all values in the input are true.
 *
 * Overload for a tuple input.
 *
 * @tparam Type of container
 * @param x Nested container of boolean inputs
 * @return Boolean indicating whether all elements are true
 */
template <typename... Types>
inline bool all(const std::tuple<Types...>& x) {
  return math::apply([](const auto&... args) {
    return all(std::vector<bool>{ all(args)... });
  }, x);
}

}  // namespace math
}  // namespace stan

#endif
