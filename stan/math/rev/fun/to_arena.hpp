#ifndef STAN_MATH_REF_FUN_TO_ARENA_HPP
#define STAN_MATH_REF_FUN_TO_ARENA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>
#include <cstring>

namespace stan {
namespace math {

/**
 * Converts given argument into a type that has any dynamic allocation on AD
 * stack.
 *
 * For types that already have this property (including scalars and
 * `var_value`s) this is a no-op.
 * @tparam T type of scalar
 * @param a argument
 * @return argument
 */
template <typename T, require_same_t<T, arena_t<T>>* = nullptr>
inline T to_arena(const T& a) {
  // intentionally making a copy (not using forwarding or returning references)
  // as these types are cheap to copy and any object referenced by an input
  // reference might go out of scope before the returned value is used
  return a;
}

/**
 * Converts given argument into a type that has any dynamic allocation on AD
 * stack.
 *
 * For std vectors that have data already on AD stack this is a shallow copy.
 * @tparam T type of scalar
 * @param a argument
 * @return argument
 */
template <typename T>
inline std::vector<T, arena_allocator<T>> to_arena(
    const std::vector<T, arena_allocator<T>>& a) {
  // What we want to do here is the same as moving input into output, except
  // that we want input to be left unchanged. With any normal allocator that
  // lead to deallocating memory twice (probably segfaulting). However,
  // dealocation with `AD_stack_allocator` is a no-op, so we can do that.
  std::vector<T, arena_allocator<T>> res;
  std::memcpy(static_cast<void*>(&res), static_cast<const void*>(&a),
              sizeof(std::vector<T, arena_allocator<T>>));
  return res;
}

/**
 * Converts given argument into a type that has any dynamic allocation on AD
 * stack.
 *
 * Converts eigen types to `arena_matrix`.
 * @tparam T type of argument
 * @param a argument
 * @return argument copied/evaluated on AD stack
 */
template <typename T, require_eigen_t<T>* = nullptr,
          require_not_same_t<T, arena_t<T>>* = nullptr>
inline arena_t<T> to_arena(const T& a) {
  return {a};
}

/**
 * Converts given argument into a type that has any dynamic allocation on AD
 * stack.
 *
 * Std vectors are copied into another std vector with custom allocator that
 * uses AD stack.
 *
 * This overload works on vectors with simple scalars that don't need to be
 * converthed themselves.
 *
 * @tparam T type of argument
 * @param a argument
 * @return argument copied on AD stack
 */
template <typename T, require_same_t<T, arena_t<T>>* = nullptr>
inline arena_t<std::vector<T>> to_arena(const std::vector<T>& a) {
  return {a.begin(), a.end()};
}

/**
 * Converts given argument into a type that has any dynamic allocation on AD
 * stack.
 *
 * Std vectors are copied into another std vector with custom allocator that
 * uses AD stack.
 *
 * This overload works on vectors with scalars that also need conversion.
 *
 * @tparam T type of argument
 * @param a argument
 * @return argument copied on AD stack
 */
template <typename T, require_not_same_t<T, arena_t<T>>* = nullptr>
inline arena_t<std::vector<T>> to_arena(const std::vector<T>& a) {
  arena_t<std::vector<T>> res;
  res.reserve(a.size());
  for (const T& i : a) {
    res.push_back(to_arena(i));
  }
  return res;
}

/**
 * If the condition is true, converts given argument into a type that has any
 * dynamic allocation on AD stack. Otherwise this is a no-op.
 *
 * @tparam T type of argument
 * @param a argument
 * @return argument copied/evaluated on AD stack
 */
template <bool Condition, typename T, std::enable_if_t<!Condition>* = nullptr>
inline T to_arena_if(T&& a) {
  return std::forward<T>(a);
}

template <bool Condition, typename T, std::enable_if_t<Condition>* = nullptr>
inline arena_t<T> to_arena_if(const T& a) {
  return to_arena(a);
}

}  // namespace math
}  // namespace stan

#endif
