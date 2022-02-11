#ifndef STAN_MATH_REV_CORE_GRAD_HPP
#define STAN_MATH_REV_CORE_GRAD_HPP

#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>
#include <stan/math/rev/core/empty_nested.hpp>
#include <stan/math/rev/core/nested_size.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Compute the gradient for all variables starting from the end of the AD tape.
 * This function does not recover memory.  The chain
 * rule is applied working down the stack from the last vari created on the
 * AD tape and then calling each vari's `chain()` method in turn.
 *
 * <p>This function computes a nested gradient only going back as far
 * as the last nesting.
 *
 * <p>This function does not recover any memory from the computation.
 *
 */
static void grad() {
  size_t end = ChainableStack::instance_->var_stack_.size();
  size_t beginning = empty_nested() ? 0 : end - nested_size();
  for (size_t i = end; i-- > beginning;) {
    ChainableStack::instance_->var_stack_[i]->chain();
  }
}

/**
 * Compute the gradient for all variables starting from the
 * specified root variable implementation.  Does not recover
 * memory.  This chainable variable's adjoint is initialized using
 * the method <code>init_dependent()</code> and then the chain
 * rule is applied working down the stack from this vari and
 * calling each vari's <code>chain()</code> method in turn.
 *
 * <p>This function computes a nested gradient only going back as far
 * as the last nesting.
 *
 * <p>This function does not recover any memory from the computation.
 *
 * @param vi Variable implementation for root of partial
 * derivative propagation.
 */
template <typename Vari>
static void grad(Vari* vi) {
  vi->init_dependent();
  grad();
}

}  // namespace math
}  // namespace stan

#endif
