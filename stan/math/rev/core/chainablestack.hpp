#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackStorage<vari, chainable_alloc> ChainableStack;

ChainableStack& chainable_stack() {
  static ChainableStack stack;
  return stack;
}

}  // namespace math
}  // namespace stan
#endif
