#ifndef STAN_MATH_REV_CORE_ARENA_ALLOCATOR_HPP
#define STAN_MATH_REV_CORE_ARENA_ALLOCATOR_HPP

#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {

/**
 * std library compatible allocator that uses AD stack.
 * @tparam T type of scalar
 */
template <typename T>
struct arena_allocator {
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using size_type = size_t;
  using difference_type = ptrdiff_t;
  using propagate_on_container_move_assignment = std::true_type;

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
   * Returns the maximum theoretically possible value of n,
   *  for which the call allocate(n, 0) could succeed.
   */
  static constexpr size_type max_size() noexcept {
    return std::numeric_limits<size_type>::max() / sizeof(T);
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
