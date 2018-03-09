#ifndef STAN_MATH_REV_CORE_ALLOCATOR_HPP
#define STAN_MATH_REV_CORE_ALLOCATOR_HPP

#include <stan/math/rev/core/chainablestack.hpp>

template <class T>
class autodiff_allocator {
 public:
  using value_type = T;

  inline autodiff_allocator() noexcept {}
  template <class U>
  autodiff_allocator(autodiff_allocator<U> const&) noexcept {}

  value_type*  // Use pointer if pointer is not a value_type*
  allocate(std::size_t n) {
    using stan::math::ChainableStack;
    return static_cast<value_type*>(
        ChainableStack::memalloc_.alloc(n * sizeof(value_type)));
  }

  void deallocate(value_type* p, std::size_t) noexcept {}
};

template <class T, class U>
bool operator==(autodiff_allocator<T> const&,
                autodiff_allocator<U> const&) noexcept {
  return true;
}

template <class T, class U>
bool operator!=(autodiff_allocator<T> const& x,
                autodiff_allocator<U> const& y) noexcept {
  return !(x == y);
}
#endif
