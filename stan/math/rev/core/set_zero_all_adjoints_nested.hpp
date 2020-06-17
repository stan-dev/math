#ifndef STAN_MATH_REV_CORE_SET_ZERO_ALL_ADJOINTS_NESTED_HPP
#define STAN_MATH_REV_CORE_SET_ZERO_ALL_ADJOINTS_NESTED_HPP

#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/empty_nested.hpp>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Reset all adjoint values in the top nested portion of the stack
 * to zero.
 *
 * It is preferred to use the <code>nested_rev_autodiff</code> class for
 * nested autodiff class as it handles recovery of memory automatically.
 */
static EIGEN_STRONG_INLINE void set_zero_all_adjoints_nested() {
  if (empty_nested()) {
    throw std::logic_error(
        "empty_nested() must be false before calling"
        " set_zero_all_adjoints_nested()");
  }
  size_t start1 = ChainableStack::instance_->nested_var_stack_sizes_.back();
  for_each_tuple([&start1](auto& x, auto& x_size) {
    const auto stack_size = x.size();
    if (stack_size > 0) {
      for (size_t i = (start1 == 0U) ? 0U : (start1 - 1); i < stack_size; ++i) {
        x[i]->set_zero_adjoint();
      }
    }
  }, ChainableStack::instance_->var_zeroing_stacks_, ChainableStack::instance_->nested_var_zeroing_stack_sizes_);
  // avoid wrap with unsigned when start1 == 0
}

}  // namespace math
}  // namespace stan
#endif
