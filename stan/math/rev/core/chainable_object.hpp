#ifndef STAN_MATH_REV_CORE_CHAINABLE_OBJECT_HPP
#define STAN_MATH_REV_CORE_CHAINABLE_OBJECT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * `chainable_object` hold another object is useful for connecting
 * the lifetime of a specific object to the chainable stack
 *
 * `chainable_object` objects should only be allocated with `new`.
 * `chainable_object` objects allocated on the stack will result
 * in a double free (`obj_` will get destructed once when the
 * chainable_object leaves scope and once when the chainable
 * stack memory is recovered).
 *
 * @tparam T type of object to hold
 */
template <typename T>
class chainable_object : public chainable_alloc {
 private:
  plain_type_t<T> obj_;

 public:
  /**
   * Construct chainable object from another object
   *
   * @tparam S type of object to hold (must have the same plain type as `T`)
   */
  template <typename S,
            require_same_t<plain_type_t<T>, plain_type_t<S>>* = nullptr>
  explicit chainable_object(S&& obj) : obj_(std::forward<S>(obj)) {}

  /**
   * Return a reference to the underlying object
   *
   * @return reference to underlying object
   */
  inline auto& get() noexcept { return obj_; }
  inline const auto& get() const noexcept { return obj_; }
};

/**
 * Store the given object in a `chainable_object` so it is destructed
 * only when the chainable stack memory is recovered and return
 * a pointer to the underlying object
 *
 * @tparam T type of object to hold
 * @param obj object to hold
 * @return pointer to object held in `chainable_object`
 */
template <typename T>
auto make_chainable_ptr(T&& obj) {
  auto ptr = new chainable_object<T>(std::forward<T>(obj));
  return &ptr->get();
}

/**
 * `unsafe_chainable_object` hold another object and is useful for connecting
 * the lifetime of a specific object to the chainable stack. This class
 * differs from `chainable_object` in that this class does not evaluate
 * expressions.
 *
 * `unsafe_chainable_object` objects should only be allocated with `new`.
 * `unsafe_chainable_object` objects allocated on the stack will result
 * in a double free (`obj_` will get destructed once when the
 * `unsafe_chainable_object` leaves scope and once when the chainable
 * stack memory is recovered).
 *
 * @tparam T type of object to hold
 */
template <typename T>
class unsafe_chainable_object : public chainable_alloc {
 private:
  std::decay_t<T> obj_;

 public:
  /**
   * Construct chainable object from another object
   *
   * @tparam S type of object to hold (must have the same plain type as `T`)
   */
  template <typename S,
            require_same_t<plain_type_t<T>, plain_type_t<S>>* = nullptr>
  explicit unsafe_chainable_object(S&& obj) : obj_(std::forward<S>(obj)) {}

  /**
   * Return a reference to the underlying object
   *
   * @return reference to underlying object
   */
  inline auto& get() noexcept { return obj_; }
  inline const auto& get() const noexcept { return obj_; }
};

/**
 * Store the given object in a `chainable_object` so it is destructed
 * only when the chainable stack memory is recovered and return
 * a pointer to the underlying object This function
 * differs from `make_chainable_object` in that this class does not evaluate
 * expressions.
 *
 * @tparam T type of object to hold
 * @param obj object to hold
 * @return pointer to object held in `chainable_object`
 */
template <typename T>
auto make_unsafe_chainable_ptr(T&& obj) {
  auto ptr = new unsafe_chainable_object<T>(std::forward<T>(obj));
  return &ptr->get();
}

}  // namespace math
}  // namespace stan
#endif
