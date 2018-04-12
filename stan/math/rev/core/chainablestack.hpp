#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackStorage<vari, chainable_alloc> ChainableStack;

#ifdef STAN_THREADS
thread_local
#endif
static ChainableStack chainable_stack;

}  // namespace math
}  // namespace stan
#endif
