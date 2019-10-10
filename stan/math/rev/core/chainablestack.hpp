#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

using ChainableStack = AutodiffStackSingleton<vari, chainable_alloc>;

// Helper struct to access the underlying stack allocator.
struct stack_mem {
  template <typename T>
  inline static T* alloc_array(size_t n) {
    return ChainableStack::instance_->memalloc_.alloc_array<T>(n);
  }
  inline static void* alloc(size_t len)  {
    return ChainableStack::instance_->memalloc_.alloc(len);
  }

  inline static void recover_all()  {
    return ChainableStack::instance_->memalloc_.recover_all();
  }
  inline static void start_nested()  {
    return ChainableStack::instance_->memalloc_.start_nested();
  }
  inline static void recover_nested()  {
    return ChainableStack::instance_->memalloc_.recover_nested();
  }
  inline static void free_all()  {
    return ChainableStack::instance_->memalloc_.free_all();
  }
  inline static size_t bytes_allocated()  {
    return ChainableStack::instance_->memalloc_.bytes_allocated();
  }
  inline static bool in_stack(const void* ptr)  {
    return ChainableStack::instance_->memalloc_.in_stack(ptr);
  }
};


}  // namespace math
}  // namespace stan
#endif
