#ifndef STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_TBB_HPP
#define STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_TBB_HPP

#include <stan/math/parallel/get_num_threads.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/task_scheduler_observer.h>

#include <mutex>
#include <iostream>

namespace stan {
namespace math {
namespace {

static std::mutex cout_mutex;

// initialize AD tape & pin thread to hardware thread
class ad_tape_observer : public tbb::task_scheduler_observer {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;

 public:
  // affinity_mask_t m_mask; // HW affinity mask to be used with an arena

  ad_tape_observer() : tbb::task_scheduler_observer() {
    std::cout << "Observer constructing..." << std::endl;
    observe(true);  // activate the observer
    std::cout << "Observer created" << std::endl;
  }

  ~ad_tape_observer() { std::cout << "Observer destructed" << std::endl; }

  /*override*/ void on_scheduler_entry(bool worker) {
    stan::math::ChainableStack::init();
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    if (worker)
      std::cout << "a worker thread is joining to do some work." << std::endl;
    else
      std::cout << "the main thread is joining to do some work." << std::endl;
  }
  /*override*/ void on_scheduler_exit(bool worker) {
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    chainablestack_t* ad_tape = stan::math::ChainableStack::instance_;
    delete ad_tape;
    stan::math::ChainableStack::instance_ = nullptr;
    if (worker)
      std::cout << "a worker thread is finished to do some work." << std::endl;
    else
      std::cout << "the main thread is finished to do some work." << std::endl;
  }
};

// the tbb_ressources is only there to ensure that the destruction is
// in correct order (first the observer, then the scheduler)
struct tbb_ressources {
  tbb_ressources()
      : tbb_scheduler(stan::math::internal::get_num_threads()),
        tape_initializer() {}

  tbb::task_scheduler_init tbb_scheduler;
  ad_tape_observer tape_initializer;
};

static tbb_ressources tbb_init;

// const ChainableStack::AutodiffStackStorage* __chainable_stack
//    = ChainableStack::init();

}  // namespace
}  // namespace math
}  // namespace stan
#endif
