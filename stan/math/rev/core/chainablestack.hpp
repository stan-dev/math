#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>

#include <iostream>

namespace stan {
namespace math {

class vari;
class chainable_alloc;

typedef AutodiffStackSingleton<vari, chainable_alloc> ChainableStack;

/**
 * Instantiates an instance of the ChainableStack if not already
 * initialized. This function must be called before any autodiff
 * variable get's instantiated within any thread which performs
 * reverse mode autodiff operations.
 */
/*
static inline ChainableStack::AutodiffStackStorage* init() {
#ifdef STAN_THREADS
  return ChainableStack::init();
  if (ChainableStack::instance_ == nullptr)
    ChainableStack::instance_ = new ChainableStack::AutodiffStackStorage();
#endif
  return ChainableStack::instance_;
}
*/

namespace {

struct ad_tape_initializer {
  ad_tape_initializer() {
    // std::cout << "Initializer created" << std::endl;
    tape_ = stan::math::ChainableStack::init();
  }
  ~ad_tape_initializer() {
    // std::lock_guard<std::mutex> cout_lock(cout_mutex_init);
    // std::cout << "Initializer destructed" << std::endl;
  }

  typedef ChainableStack::AutodiffStackStorage* tape_ptr_t;
  thread_local static tape_ptr_t tape_;
};

thread_local ad_tape_initializer::tape_ptr_t ad_tape_initializer::tape_
    = ChainableStack::init();
ad_tape_initializer global_initializer;

}  // namespace

}  // namespace math
}  // namespace stan
#endif
