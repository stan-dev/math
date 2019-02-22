#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
#include <iostream>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackSingleton<vari, chainable_alloc> ChainableStack;

// namespace {

/**/
struct ChainableStackInit {
  ChainableStackInit() {
    std::cout << "ChainableStackInit instance." << std::endl;
    ChainableStack::init_instance();
  }
  ~ChainableStackInit() {}

  typedef ChainableStack::AutodiffStackStorage* chainable_stack_ptr_t;

  thread_local static chainable_stack_ptr_t instance_;
};

#ifdef STAN_THREADS
thread_local
#endif
    ChainableStackInit::chainable_stack_ptr_t ChainableStackInit::instance_
    = ChainableStack::init_instance();

static ChainableStackInit thread_initializer;
/**/
//}  // namespace

}  // namespace math
}  // namespace stan
#endif
