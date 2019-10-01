#ifndef STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_TBB_HPP
#define STAN_MATH_REV_CORE_INIT_CHAINABLESTACK_TBB_HPP

// TODO(SW): remove STAN_THREADS guard once Intel TBB is fully
// mandatory
#ifdef STAN_THREADS

#include <stan/math/rev/core/chainablestack.hpp>

#include <tbb/task_scheduler_observer.h>

#include <mutex>
#include <iostream>
#include <unordered_map>
#include <thread>
#include <utility>
#include <tuple>

namespace stan {
namespace math {

/**
 * TBB observer object which is a callback hook called whenever a new
 * thread enters the threadpool managed by the TBB. This hook ensures
 * that each worker thread has an initialized AD tape ready for use.
 */
class ad_tape_observer : public tbb::task_scheduler_observer {
  using ad_map = std::unordered_map<std::thread::id, ChainableStack*>;

 public:
  ad_tape_observer() : tbb::task_scheduler_observer(), thread_tape_map_() {
    observe(true);  // activates the observer
  }

  ~ad_tape_observer() {
    for (auto& elem : thread_tape_map_) {
      delete elem.second;
    }
    thread_tape_map_.clear();
  }

  void on_scheduler_entry(bool worker) {
    const std::thread::id thread_id = std::this_thread::get_id();
    if (thread_tape_map_.find(thread_id) == thread_tape_map_.end()) {
      ad_map::iterator insert_elem;
      bool status = false;
      std::tie(insert_elem, status)
          = thread_tape_map_.emplace(ad_map::value_type{thread_id, nullptr});
      insert_elem->second = new ChainableStack();
    }
  }

  void on_scheduler_exit(bool worker) {
    auto elem = thread_tape_map_.find(std::this_thread::get_id());
    if (elem != thread_tape_map_.end()) {
      delete elem->second;
      thread_tape_map_.erase(elem);
    }
  }

 private:
  ad_map thread_tape_map_;
};

namespace {

static ad_tape_observer global_observer;

}  // namespace
}  // namespace math
}  // namespace stan

#endif  // STAN_THREADS

#endif
