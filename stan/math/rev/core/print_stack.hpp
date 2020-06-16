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
}

}  // namespace math
}  // namespace stan
#endif
