#ifndef STAN_MATH_REV_CORE_ARENA_ALLOCATOR_HPP
#define STAN_MATH_REV_CORE_ARENA_ALLOCATOR_HPP

#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {

/**
 * std library compatible allocator that uses AD stack.
 * @tparam T type of scalar
 *
 * @warning The type T needs to be either trivially destructible or the dynamic
allocations needs to be managed by the arena_allocator.
 * For example this works: @code{.cpp}
using my_matrix = std::vector<std::vector<double,
stan::math::arena_allocator<double>>,
stan::math::arena_allocator<std::vector<double,
stan::math::arena_allocator<double>>>>;@endcode
 *
 */
template <typename T>
struct arena_allocator {
  using value_type = T;

  arena_allocator() = default;

  arena_allocator(const arena_allocator& rhs) = default;

  template <class U>
  arena_allocator(const arena_allocator<U>& rhs) {}

  /**
   * Allocates space for `n` items of type `T`.
   *
   * @param n number of items to allocate space for
   * @return pointer to allocated space
   */
  T* allocate(std::size_t n) {
    return ChainableStack::instance_->memalloc_.alloc_array<T>(n);
  }

  /**
   * No-op. Memory is deallocated by calling `recover_memory()`.
   */
  void deallocate(T* /*p*/, std::size_t /*n*/) noexcept {}

  /**
   * Equality comparison operator.
   * @return true
   */
  constexpr bool operator==(const arena_allocator&) const noexcept {
    return true;
  }
  /**
   * Inequality comparison operator.
   * @return false
   */
  constexpr bool operator!=(const arena_allocator&) const noexcept {
    return false;
  }
};

}  // namespace math
}  // namespace stan

#endif
