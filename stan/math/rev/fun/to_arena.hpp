#ifndef STAN_MATH_REV_FUN_TO_ARENA_HPP
#define STAN_MATH_REV_FUN_TO_ARENA_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <vector>
#include <cstring>

namespace stan {
namespace math {

/**
 * Converts given argument into a type that either has any dynamic allocation on
 * AD stack or schedules its destructor to be called when AD stack memory is
 * recovered.
 *
 * This overload is used for kernel generator expressions. It also handles any
 * other types that do not have a special overload for them.
 *
 * @tparam Exprcl type of scalar
 * @param a argument
 * @return argument
 */
template <typename Exprcl,
          require_not_same_t<Exprcl, arena_t<Exprcl>>* = nullptr,
          require_not_container_t<Exprcl>* = nullptr,
          require_not_matrix_cl_t<Exprcl>* = nullptr>
inline arena_t<Exprcl> to_arena(Exprcl&& a) {
  return std::forward<Exprcl>(a);
}

/**
 * Converts given argument into a type that either has any dynamic allocation on
 * AD stack or schedules its destructor to be called when AD stack memory is
 * recovered.
 *
 * For types that already have this property (including scalars and
 * \ref stan::math::var_value ) this is a no-op.
 *
 * Passing in a lvalue reference to objects not using AD stack, such as a
 * `matrix_cl` is inefficient as they need to be copied in this case.
 * @tparam ArenaT type of scalar
 * @param a argument
 * @return argument
 */
template <typename ArenaT, require_same_t<ArenaT, arena_t<ArenaT>>* = nullptr,
          require_not_matrix_cl_t<ArenaT>* = nullptr,
          require_not_std_vector_t<ArenaT>* = nullptr>
inline std::remove_reference_t<ArenaT> to_arena(ArenaT&& a) {
  // intentionally never returning a reference. If an object is just
  // referenced it will likely go out of scope before it is used.
  return std::forward<ArenaT>(a);
}

/**
 * Converts given argument into a type that either has any dynamic allocation on
 * AD stack or schedules its destructor to be called when AD stack memory is
 * recovered.
 *
 * Converts eigen types to \ref stan::math::arena_matrix.
 * @tparam Eig type of argument
 * @param a argument
 * @return argument copied/evaluated on AD stack
 */
template <typename Eig, require_eigen_t<Eig>* = nullptr,
          require_not_same_t<Eig, arena_t<Eig>>* = nullptr>
inline arena_t<Eig> to_arena(const Eig& a) {
  return arena_t<Eig>(a);
}

/**
 * Converts given argument into a type that either has any dynamic allocation on
 * AD stack or schedules its destructor to be called when AD stack memory is
 * recovered.
 *
 * For std vectors that have data already on AD stack this is a shallow copy.
 * @tparam Scalar type of scalar
 * @param a argument
 * @return argument
 */
template <typename Scalar>
inline std::vector<Scalar, arena_allocator<Scalar>> to_arena(
    const std::vector<Scalar, arena_allocator<Scalar>>& a) {
  // What we want to do here is the same as moving input into output, except
  // that we want input to be left unchanged. With any normal allocator that
  // lead to deallocating memory twice (probably segfaulting). However,
  // dealocation with `arena_allocator` is a no-op, so we can do that.
  std::vector<Scalar, arena_allocator<Scalar>> res;
  std::memcpy(static_cast<void*>(&res), static_cast<const void*>(&a),
              sizeof(std::vector<Scalar, arena_allocator<Scalar>>));
  return res;
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
 * @tparam ArenaT type of argument
 * @param a argument
 * @return argument copied on AD stack
 */
template <typename ArenaT, require_same_t<ArenaT, arena_t<ArenaT>>* = nullptr>
inline arena_t<std::vector<ArenaT>> to_arena(const std::vector<ArenaT>& a) {
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
 * dynamic allocation on AD stack. Otherwise returns the argument
 *
 * @tparam T1 type of argument
 * @param a argument
 * @return argument copied/evaluated on AD stack
 */
template <bool Condition, typename T1, std::enable_if_t<!Condition>* = nullptr>
inline T1 to_arena_if(T1&& a) {
  return std::forward<T1>(a);
}

template <bool Condition, typename T2, std::enable_if_t<Condition>* = nullptr>
inline arena_t<T2> to_arena_if(const T2& a) {
  return to_arena(a);
}

}  // namespace math
}  // namespace stan

#endif
