#ifndef STAN_MATH_REV_CORE_PRINT_STACK_HPP
#define STAN_MATH_REV_CORE_PRINT_STACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <ostream>

namespace stan {
namespace math {

/**
 * Prints the autodiff variable stack. This function
 * is used for debugging purposes.
 *
 * Only works if all members of stack are vari* as it
 * casts to vari*.
 *
 * @param o ostream to modify
 */
inline void print_stack(std::ostream& o) {
  o << "STACK, size=" << ChainableStack::instance_->var_stack_.size()
    << std::endl;
  // TODO(carpenter): this shouldn't need to be cast any more
  for (size_t i = 0; i < ChainableStack::instance_->var_stack_.size(); ++i) {
    boost::variant2::visit(
        [&o, &i](auto x) {
          o << i << "  " << x << "  " << x->val_ << " : " << x->adj_
            << std::endl;
        },
        ChainableStack::instance_->var_stack_[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif
