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
   * Constructs an object of type T in allocated uninitialized storage pointed
   * to by p, using placement-new
   * @param p a pointer to type `T`
   * @param val The value to use during construction.
   */
  static inline void construct(pointer p, const_reference val) {
    new (static_cast<void*>(p)) T(val);
  }
  /**
   * Constructs an object of type U in allocated uninitialized storage pointed
   * to by p, using placement-new
   * @tparam U The type to construct.
   * @tparam Args Parameter pack of types used in constructor.
   * @param p a pointer of type `U`
   * @param args Parameter pack of objects used in constructor.
   */
  template <class U, class... Args>
  void construct(U* p, Args&&... args) {
    new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
  }

  /**
   * Use the allocator `a` to construct an object of type `U` in allocated
   *  uninitialized storage pointed to by `p`, using placement-new.
   *  This function is used by the standard library containers when inserting,
   *  copying, or moving elements.
   * @tparam U The type to construct.
   * @tparam Args Parameter pack of types used in constructor.
   * @param a An arena allocator that performs the construction.
   * @param p a pointer of type `U`
   * @param args Parameter pack of objects used in constructor.
   */
  template <class U, class... Args>
  static void construct(arena_allocator<U>& a, U* p, Args&&... args) {
    a.construct(p, args...);
  }

  /**
   * No-op. If anything needs destroyed on the stack this happens
   *  when calling `recover_memory()`.
   */
  static inline void destroy(pointer p) {}

  /**
   * Returns the actual address of x even in presence of overloaded operator&
   */
  inline pointer address(reference x) const noexcept { return &x; }

  /**
   * Returns the actual address of x even in presence of overloaded operator&
   */
  inline const_pointer address(const_reference x) const noexcept { return &x; }

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
