#ifndef STAN_MATH_REV_CORE_SET_ZERO_ALL_ADJOINTS_HPP
#define STAN_MATH_REV_CORE_SET_ZERO_ALL_ADJOINTS_HPP

#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {

/**
 * Reset all adjoint values in the stack to zero.
 */
static EIGEN_STRONG_INLINE void set_zero_all_adjoints() {
  for_each_tuple([](auto& x) {
    for (auto& xx : x) {
      xx->set_zero_adjoint();
    }
    return 0;
  }, ChainableStack::instance_->var_zeroing_stacks_);
}

}  // namespace math
}  // namespace stan
#endif
