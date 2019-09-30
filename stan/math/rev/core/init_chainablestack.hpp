#ifndef STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_HPP

#include <stan/math/rev/core/chainablestack.hpp>

#include <tbb/task_scheduler_observer.h>

#include <mutex>
#include <iostream>
#include <unordered_map>
#include <thread>
#include <utility>

namespace stan {
namespace math {

static std::mutex cout_mutex;

// initialize AD tape & pin thread to hardware thread
class ad_tape_observer : public tbb::task_scheduler_observer {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;

 public:

  ad_tape_observer() : tbb::task_scheduler_observer(), ad_tape_() {
    std::cout << "Observer constructing..." << std::endl;
    observe(true);  // activate the observer
    std::cout << "Observer created" << std::endl;
  }

  ~ad_tape_observer() { std::cout << "Observer destructed" << std::endl; }

  /*override*/ void on_scheduler_entry(bool worker) {
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    const std::thread::id thread_id = std::this_thread::get_id();
    if(ad_tape_.find(thread_id) != ad_tape_.end()) {
      std::cout << "Initialising a thread AD tape." << std::endl;
      ad_tape_.emplace(std::pair<const std::thread::id,ChainableStack>{std::this_thread::get_id(), ChainableStack()});
    }
    if (worker)
      std::cout << "a worker thread is joining to do some work." << std::endl;
    else
      std::cout << "the main thread is joining to do some work." << std::endl;
  }
  /*override*/ void on_scheduler_exit(bool worker) {
    std::lock_guard<std::mutex> cout_lock(cout_mutex);
    const std::thread::id thread_id = std::this_thread::get_id();
    ad_tape_.erase(thread_id);
    if (worker)
      std::cout << "a worker thread is finished to do some work." << std::endl;
    else
      std::cout << "the main thread is finished to do some work." << std::endl;
  }

 private:
  std::unordered_map<std::thread::id,ChainableStack> ad_tape_;
};

namespace {

/**
 * Initializes the AD stack for the main process. See
 * autodiffstackstorage.hpp for more explanations.
 */
//ChainableStack global_stack_instance_init;
static ad_tape_observer global_observer ;

}  // namespace
}  // namespace math
}  // namespace stan
#endif
