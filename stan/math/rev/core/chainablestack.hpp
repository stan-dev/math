#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackStorage<vari, chainable_alloc> ChainableStack_t;

#ifndef STAN_THREADS
static ChainableStack_t ChainableStack;
#else
static thread_local ChainableStack_t ChainableStack;
#endif

}  // namespace math
}  // namespace stan
#endif
