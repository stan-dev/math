#ifndef STAN_MATH_PRIM_SCAL_INIT_THREADPOOL_TBB_HPP
#define STAN_MATH_PRIM_SCAL_INIT_THREADPOOL_TBB_HPP

// TODO(SW): remove STAN_THREADS guard once Intel TBB is fully
// mandatory
#ifdef STAN_THREADS

#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <boost/lexical_cast.hpp>

#include <tbb/task_scheduler_init.h>

#include <cstdlib>

namespace stan {
namespace math {

/**
 * Initialize the Intel TBB threadpool through the
 * tbb::task_scheduler_init object. In case an instance of the
 * tbb::task_scheduler_object is instantiated, then this has no effect
 * which is the default behavior of the Intel TBB library.
 *
 * The maximal number of threads is read from the environment variable
 * STAN_NUM_THREADS if use_env = true, which is the default. If
 * use_env = false, then the max_threads argument is passed into the
 * constructor of tbb::task_scheduler_init as the max_threads
 * argument. The default max_threads corresponds to the default of the
 * Intel TBB.
 *
 * The STAN_NUM_THREADS variable follows these conventions:
 *
 * - STAN_NUM_THREADS is not defined => max_threads=1
 * - STAN_NUM_THREADS is positive => max_threads is set to the
 *   specified number
 * - STAN_NUM_THREADS is set to -1 => max_threads is set by the Intel
 *   TBB to the maximal hardware concurency
 * - STAN_NUM_THREADS < -1, STAN_NUM_THREADS = 0 or STAN_NUM_THREADS is
 *   not numeric => throws an exception
 *
 * The function returns as a boolean if the tbb::task_scheduler_init
 * has become active once constructed. The instance may not become
 * active if there is already another instance created before calling
 * this function.
 *
 * @param use_env determines if max_threads is determined from the
 * environment variable STAN_NUM_THREADS; defaults to true
 * @param max_threads maximal number of threads argument of the
 * tbb::task_scheduler_init constructor (overridden if use_env=true);
 * defaults to automatic selection of number of threads by the TBB
 * @param stack_size sets the stack size of each thread; the default 0
 * let's the TBB choose the stack size
 * @return active status of static tbb::task_scheduler_init
 * @throws std::runtime_error if the value of STAN_NUM_THREADS env. variable
 * is invalid
 */
inline bool init_threadpool_tbb(bool use_env = true,
                                int max_threads
                                = tbb::task_scheduler_init::automatic,
                                tbb::stack_size_type stack_size = 0) {
  int tbb_max_threads = max_threads;
  if (use_env) {
    const char* env_stan_num_threads = std::getenv("STAN_NUM_THREADS");
    if (env_stan_num_threads != nullptr) {
      try {
        const int env_max_threads
            = boost::lexical_cast<int>(env_stan_num_threads);
        if (env_max_threads > 0)
          tbb_max_threads = env_max_threads;
        else if (env_max_threads == -1)
          tbb_max_threads = tbb::task_scheduler_init::automatic;
        else
          invalid_argument("init_threadpool_tbb", "STAN_NUM_THREADS",
                           env_stan_num_threads,
                           "The STAN_NUM_THREADS environment variable is '",
                           "' but it must be positive or -1");
      } catch (boost::bad_lexical_cast) {
        invalid_argument("init_threadpool_tbb", "STAN_NUM_THREADS",
                         env_stan_num_threads,
                         "The STAN_NUM_THREADS environment variable is '",
                         "' but it must be a positive number or -1");
      }
    }
  }

  static tbb::task_scheduler_init tbb_scheduler(tbb_max_threads, stack_size);

  return tbb_scheduler.is_active();
}

}  // namespace math
}  // namespace stan

#endif

#endif
