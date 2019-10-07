#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

using ChainableStack = AutodiffStackSingleton<vari, chainable_alloc>;

// Helper struct to access the underlying stack allocator.
struct stack_mem_adapter {
  template <typename T>
  inline T* alloc_array(size_t n) const {
    return ChainableStack::instance_->memalloc_.alloc_array<T>(n);
  }
  inline void* alloc(size_t len) const {
    return ChainableStack::instance_->memalloc_.alloc<double>(len);
  }

  inline void recover_all() const {
    return ChainableStack::instance_->memalloc_.recover_all();
  }
  inline void start_nested() const {
    return ChainableStack::instance_->memalloc_.start_nested();
  }
  inline void recover_nested() const {
    return ChainableStack::instance_->memalloc_.recover_nested();
  }
  inline void free_all() const {
    return ChainableStack::instance_->memalloc_.free_all();
  }
  inline size_t bytes_allocated() const {
    return ChainableStack::instance_->memalloc_.bytes_allocated();
  }
  inline bool in_stack(const void* ptr) const {
    return ChainableStack::instance_->memalloc_.in_stack(ptr);
  }
};

constexpr stack_mem_adapter stack_mem;

}  // namespace math
}  // namespace stan
#endif
