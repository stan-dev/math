#ifndef STAN_MATH_PRIM_FUN_ANY_HPP
#define STAN_MATH_PRIM_FUN_ANY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Return true if any values in the input are true.
 *
 * Overload for a single boolean input
 *
 * @tparam T Any type convertible to `bool`
 * @param x boolean input
 * @return The input unchanged
 */
template <typename T, require_t<std::is_convertible<T, bool>>* = nullptr>
constexpr inline bool any(T x) {
  return x;
}

/**
 * Return true if any values in the input are true.
 *
 * Overload for Eigen types
 *
 * @tparam ContainerT A type derived from `Eigen::EigenBase` that has an
 *                      `integral` scalar type
 * @param x Eigen object of boolean inputs
 * @return Boolean indicating whether any elements are true
 */
template <typename ContainerT,
          require_eigen_st<std::is_integral, ContainerT>* = nullptr>
inline bool any(const ContainerT& x) {
  return x.any();
}

// Forward-declaration for correct resolution of any(std::vector<std::tuple>)
template <typename... Types>
inline bool any(const std::tuple<Types...>& x);

/**
 * Return true if any values in the input are true.
 *
 * Overload for a std::vector/nested inputs. The Eigen::Map/apply_vector_unary
 * approach cannot be used as std::vector<bool> types do not have a .data()
 * member and are not always stored contiguously.
 *
 * @tparam InnerT Type within std::vector
 * @param x Nested container of boolean inputs
 * @return Boolean indicating whether any elements are true
 */
template <typename InnerT>
inline bool any(const std::vector<InnerT>& x) {
  return std::any_of(x.begin(), x.end(), [](const auto& i) { return any(i); });
}

/**
 * Return true if any values in the input are true.
 *
 * Overload for a tuple input.
 *
 * @tparam Types of items within tuple
 * @param x Tuple of boolean scalar-type elements
 * @return Boolean indicating whether any elements are true
 */
template <typename... Types>
inline bool any(const std::tuple<Types...>& x) {
  bool any_true = false;
  math::for_each(
      [&any_true](const auto& i) {
        any_true = any_true || any(i);
        return;
      },
      x);
  return any_true;
}

}  // namespace math
}  // namespace stan

#endif
