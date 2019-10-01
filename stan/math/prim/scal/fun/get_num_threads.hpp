#ifndef STAN_MATH_PRIM_SCAL_FUN_GET_NUM_THREADS_HPP
#define STAN_MATH_PRIM_SCAL_FUN_GET_NUM_THREADS_HPP

#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <boost/lexical_cast.hpp>

#include <cstdlib>
#include <thread>

namespace stan {
namespace math {
namespace internal {

/**
 * Get number of threads to use for num_jobs jobs. The function uses
 * the environment variable STAN_NUM_THREADS and follows these
 * conventions:
 *
 * - STAN_NUM_THREADS is not defined => num_threads=1
 * - STAN_NUM_THREADS is positive => num_threads is set to the
 *   specified number
 * - STAN_NUM_THREADS is set to -1 => num_threads is the number of
 *   available cores on the machine
 * - STAN_NUM_THREADS < -1, STAN_NUM_THREADS = 0 or STAN_NUM_THREADS is
 *   not numeric => throws an exception
 *
 * Should num_threads exceed the number of jobs, then num_threads will
 * be set equal to the number of jobs.
 *
 * @param num_jobs number of jobs
 * @return number of threads to use
 * @throws std::runtime_error if the value of STAN_NUM_THREADS env. variable
 * is invalid
 */
// TODO(SW):
// 1. get rid of num_jobs argument which is not needed to find
// out about number of threads we need to use.
// 2. Consider to switch from std::thread::hardware_concurrency() to
// tbb::task_scheduler_init::default_num_threads() once TBB is
// mandatory.
// 3. pull out of internal?
inline int get_num_threads(int num_jobs) {
  int num_threads = 1;
#ifdef STAN_THREADS
  const char* env_stan_num_threads = std::getenv("STAN_NUM_THREADS");
  if (env_stan_num_threads != nullptr) {
    try {
      const int env_num_threads
          = boost::lexical_cast<int>(env_stan_num_threads);
      if (env_num_threads > 0)
        num_threads = env_num_threads;
      else if (env_num_threads == -1)
        num_threads = std::thread::hardware_concurrency();
      else
        invalid_argument("get_num_threads(int)", "STAN_NUM_THREADS",
                         env_stan_num_threads,
                         "The STAN_NUM_THREADS environment variable is '",
                         "' but it must be positive or -1");
    } catch (boost::bad_lexical_cast) {
      invalid_argument("get_num_threads(int)", "STAN_NUM_THREADS",
                       env_stan_num_threads,
                       "The STAN_NUM_THREADS environment variable is '",
                       "' but it must be a positive number or -1");
    }
  }
  if (num_threads > num_jobs)
    num_threads = num_jobs;
#endif
  return num_threads;
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
