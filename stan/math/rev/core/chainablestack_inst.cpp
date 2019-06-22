#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_INST_CPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_INST_CPP

#warning "Init is part of header files! You do not need this!"

#include <stan/math/parallel/get_num_threads.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/task_scheduler_observer.h>

#include <mutex>
#include <iostream>

namespace stan {
namespace math {
namespace internal {

static tbb::task_scheduler_init task_scheduler(
    stan::math::internal::get_num_threads());

static std::mutex cout_mutex;

class ad_tape_observer : public tbb::task_scheduler_observer {
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
    if (worker)
      std::cout << "a worker thread is finished to do some work." << std::endl;
    else
      std::cout << "the main thread is finished to do some work." << std::endl;
  }
};

static ad_tape_observer tape_initializer;

// const ChainableStack::AutodiffStackStorage* __chainable_stack
//    = ChainableStack::init();

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
