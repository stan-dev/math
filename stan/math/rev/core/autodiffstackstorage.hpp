#ifndef STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP
#define STAN_MATH_REV_CORE_AUTODIFFSTACKSTORAGE_HPP

#include <stan/math/memory/stack_alloc.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename ChainableT, typename ChainableAllocT>
struct AutodiffStackStorage {

  std::vector<ChainableT*> var_stack_;
  std::vector<ChainableT*> var_nochain_stack_;
  std::vector<ChainableAllocT*> var_alloc_stack_;
  stack_alloc memalloc_;

  // nested positions
  std::vector<size_t> nested_var_stack_sizes_;
  std::vector<size_t> nested_var_nochain_stack_sizes_;
  std::vector<size_t> nested_var_alloc_stack_starts_;
};

}  // namespace math
}  // namespace stan
#endif
