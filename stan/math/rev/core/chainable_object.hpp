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
  T obj_;

 public:
  /**
   * Construct chainable object from another object
   *
   * @tparam S type of object to hold (must be the same as `T` excluding `cost`
   * and `reference` qualitifers)
   */
  template <typename S,
            require_same_t<plain_type_t<T>, plain_type_t<S>>* = nullptr>
  chainable_object(S&& obj) : obj_(std::forward<S>(obj)) {}

  /**
   * Return a reference to the underlying object
   *
   * @return reference to underlying object
   */
  inline T& get() noexcept { return obj_; }
};

/**
 * Store the given object in a `chainable_object` so it is destructed
 * only when the chainable stack memory is recovered
 *
 * @tparam T type of object to hold
 * @param object to hold
 * @return pointer to `chainable_object` holding input
 */
template <typename T>
auto make_chainable_ptr(T&& obj) {
  auto ptr = new chainable_object<std::decay_t<T>>(std::forward<T>(obj));
  return ptr;
}

}  // namespace math
}  // namespace stan
#endif
