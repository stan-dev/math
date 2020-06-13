#ifndef STAN_MATH_REV_CORE_SET_ZERO_ALL_ADJOINTS_NESTED_HPP
#define STAN_MATH_REV_CORE_SET_ZERO_ALL_ADJOINTS_NESTED_HPP

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
inline void set_zero_all_adjoints_nested() {
  if (empty_nested()) {
    throw std::logic_error(
        "empty_nested() must be false before calling"
        " set_zero_all_adjoints_nested()");
  }
  size_t start1 = ChainableStack::instance_->nested_var_stack_sizes_.back();
  const auto stack_size = ChainableStack::instance_->var_stack_.size();
  // avoid wrap with unsigned when start1 == 0
  for (size_t i = (start1 == 0U) ? 0U : (start1 - 1); i < stack_size; ++i) {
    boost::apply_visitor(vari_zero_adj(),
                         ChainableStack::instance_->var_stack_[i]);
  }

  size_t start2
      = ChainableStack::instance_->nested_var_nochain_stack_sizes_.back();
  const auto nochain_stack_size
      = ChainableStack::instance_->var_nochain_stack_.size();
  for (size_t i = (start2 == 0U) ? 0U : (start2 - 1); i < nochain_stack_size;
       ++i) {
    boost::apply_visitor(vari_zero_adj(),
                         ChainableStack::instance_->var_nochain_stack_[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif
