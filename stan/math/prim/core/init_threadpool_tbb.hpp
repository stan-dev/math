#ifndef STAN_MATH_PRIM_CORE_INIT_THREADPOOL_TBB_HPP
#define STAN_MATH_PRIM_CORE_INIT_THREADPOOL_TBB_HPP

#include <stan/math/prim/err/invalid_argument.hpp>

#include <boost/lexical_cast.hpp>

#include <tbb/task_scheduler_init.h>

#include <cstdlib>
#include <thread>
#include <iostream>

namespace stan {
namespace math {
namespace internal {

/**
 * @deprecated use <code>gp_exp_quad_cov</code>
 */
inline int get_num_threads() {
  int num_threads = 1;
  const char* env_stan_num_threads = std::getenv("STAN_NUM_THREADS");
  if (env_stan_num_threads != nullptr) {
    try {
      const int env_num_threads
          = boost::lexical_cast<int>(env_stan_num_threads);
      if (env_num_threads > 0) {
        num_threads = env_num_threads;
      } else if (env_num_threads == -1) {
        num_threads = std::thread::hardware_concurrency();
      } else {
        invalid_argument("get_num_threads(int)", "STAN_NUM_THREADS",
                         env_stan_num_threads,
                         "The STAN_NUM_THREADS environment variable is '",
                         "' but it must be positive or -1");
      }
    } catch (const boost::bad_lexical_cast&) {
      invalid_argument("get_num_threads(int)", "STAN_NUM_THREADS",
                       env_stan_num_threads,
                       "The STAN_NUM_THREADS environment variable is '",
                       "' but it must be a positive number or -1");
    }
  }
  return num_threads;
}

}  // namespace internal

/**
 * Initialize the Intel TBB threadpool and global scheduler through
 * the tbb::task_scheduler_init object. In case an instance of the
 * tbb::task_scheduler_object has been instantiated prior to calling
 * this function, then any subsequent initialization is ignored by the
 * Intel TBB.
 *
 * The maximum number of threads is set from the supplied argument.
 * The TBB scheduler will be activated by calling this function.
 *
 * The function returns a reference to the static
 * tbb::task_scheduler_init instance.
 *
 * @param max_threads maximum number of threads used; If a positive integer is
 * supplied, that determines the maximum number of threads used.
 * If -1, the maximum number of threads is set to the number of concurrent
 * threads supported. If -2 (default value), the deprecated STAN_NUM_THREADS
 * env. variable determines the maximum number of threads used. Other values
 * are invalid.
 * @param stack_size sets the stack size of each thread; the default 0
 * let's the TBB choose the stack size
 * @return reference to the static tbb::task_scheduler_init
 * @throws std::invalid_argument if the supplied value of max_threads is 
 * invalid of the deprecated value of the STAN_NUM_THREADS env. variable
 * is invalid.
 */
inline tbb::task_scheduler_init& init_threadpool_tbb(
    int max_threads = -2,
    tbb::stack_size_type stack_size = 0) {
  int tbb_max_threads = -1;
  if (max_threads == -2) {
    tbb_max_threads = internal::get_num_threads();
  }else if(max_threads == -1){
    tbb_max_threads = std::thread::hardware_concurrency();
  } else if(max_threads > 0){
    tbb_max_threads = max_threads;
  } else {
     invalid_argument("init_threadpool_tbb(int)", "max_threads",
                       max_threads,
                       "The max_threads argument is '",
                       "' but it must be a positive integer or -1");
  }
  static tbb::task_scheduler_init tbb_scheduler(tbb_max_threads, stack_size);

  return tbb_scheduler;
}

}  // namespace math
}  // namespace stan

#endif
