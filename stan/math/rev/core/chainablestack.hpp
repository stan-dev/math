#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackSingleton<vari, chainable_alloc> ChainableStack;

namespace {

struct ChainableStackInit {
  ChainableStackInit() { instance_ = ChainableStack::init_instance(); }
  ~ChainableStackInit() {}

  typedef ChainableStack::AutodiffStackStorage* chainable_stack_ptr_t;

#ifdef STAN_THREADS
  thread_local
#endif
      static chainable_stack_ptr_t instance_;
};

#ifdef STAN_THREADS
thread_local
#endif
    ChainableStackInit::chainable_stack_ptr_t ChainableStackInit::instance_
    = ChainableStack::init_instance();

ChainableStackInit global_initializer;

}  // namespace

}  // namespace math
}  // namespace stan
#endif
