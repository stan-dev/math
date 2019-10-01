#ifndef STAN_MATH_PRIM_SCAL_INIT_THREADPOOL_TBB_HPP
#define STAN_MATH_PRIM_SCAL_INIT_THREADPOOL_TBB_HPP

// TODO(SW): remove STAN_THREADS guards once Intel TBB is fully
// mandatory

#include <stan/math/prim/scal/fun/get_num_threads.hpp>

#ifdef STAN_THREADS
#include <tbb/task_scheduler_init.h>
#endif

namespace stan {
namespace math {

/**
 * Initialize the Intel TBB threadpool through the
 * tbb::task_scheduler_init object. In case an instance of the
 * tbb::task_scheduler_object is instantiated, then this has no effect,
 * which is the default behavior of the Intel TBB library.
 *
 * The maximal number of threads is read from the environment variable
 * STAN_NUM_THREADS using internal::get_num_threads. See conventions
 * of get_num_threads.
 *
 * The function returns as a boolean if the tbb::task_scheduler_init
 * has become active once constructed. The instance may not become
 * active if there is already another instance created before calling
 * this function.
 *
 * @param stack_size sets the stack size of each thread; the default 0
 * let's the TBB choose the stack size
 * @return active status of static tbb::task_scheduler_init
 * @throws std::runtime_error if the value of STAN_NUM_THREADS env. variable
 * is invalid
 */
inline bool init_threadpool_tbb(tbb::stack_size_type stack_size = 0) {
#ifdef STAN_THREADS
  int tbb_max_threads = internal::get_num_threads(
      tbb::task_scheduler_init::default_num_threads());

  static tbb::task_scheduler_init tbb_scheduler(tbb_max_threads, stack_size);

  return tbb_scheduler.is_active();
#else
  return false;
#endif
}

}  // namespace math
}  // namespace stan

#endif
