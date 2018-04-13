#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackStorage<vari, chainable_alloc> ChainableStack;

struct ADStack {
#ifdef STAN_THREADS
  thread_local
#endif
  static ChainableStack instance;
};

#ifdef STAN_THREADS
thread_local
#endif
ChainableStack ADStack::instance;

}  // namespace math
}  // namespace stan
#endif
