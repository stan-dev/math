#ifndef STAN_MATH_PRIM_CORE_INIT_THREADPOOL_TBB_HPP
#define STAN_MATH_PRIM_CORE_INIT_THREADPOOL_TBB_HPP

#include <stan/math/prim/err/invalid_argument.hpp>

#include <boost/lexical_cast.hpp>

#ifndef TBB_INTERFACE_NEW
#include <tbb/tbb_stddef.h>

#if TBB_VERSION_MAJOR >= 2020
#define TBB_INTERFACE_NEW
#endif
#endif

#ifdef TBB_INTERFACE_NEW
#include <tbb/global_control.h>
#include <tbb/task_arena.h>
#else
#include <tbb/task_scheduler_init.h>
#endif

#include <cstdlib>
#include <thread>

namespace stan {
namespace math {
namespace internal {

/**
 * Get number of threads to use. The function uses the environment
 * variable STAN_NUM_THREADS and follows these conventions:
 *
 * - STAN_NUM_THREADS is not defined => num_threads=1
 * - STAN_NUM_THREADS is positive => num_threads is set to the
 *   specified number
 * - STAN_NUM_THREADS is set to -1 => num_threads is the number of
 *   available cores on the machine
 * - STAN_NUM_THREADS < -1, STAN_NUM_THREADS = 0 or STAN_NUM_THREADS is
 *   not numeric => throws an exception
 *
 * @return number of threads to use
 * @throws std::invalid_argument if the value of STAN_NUM_THREADS env. variable
 * is invalid
 */
inline int get_num_threads() {
  int num_threads = 1;
#ifdef STAN_THREADS
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
#endif
  return num_threads;
}

}  // namespace internal

#ifdef TBB_INTERFACE_NEW
/**
 * Initialize the Intel TBB threadpool and global scheduler through
 * the tbb::task_arena object. In case an instance of the
 * tbb::task_scheduler_object has been instantiated prior to calling
 * this function, then any subsequent initialization is ignored by the
 * Intel TBB.
 *
 * The maximal number of threads is read from the environment variable
 * STAN_NUM_THREADS using internal::get_num_threads. See conventions
 * of get_num_threads. The TBB scheduler will be activated by calling
 * this function.
 *
 * The function returns a reference to the static
 * tbb::global_control instance.
 *
 * @param n_threads The maximum number of threads available to the tbb. If not
 *  set will search for the environment variable `STAN_NUM_THREADS`. A value of
 *  -1 will assume all detectable threads are available for the process.
 * @return reference to the static tbb::global_control
 * @throws std::runtime_error if n_threads (defaults to zero) is not provided
 * and the value of STAN_NUM_THREADS environment variable is invalid.
 */
inline tbb::task_arena& init_threadpool_tbb(int n_threads = 0) {
  int tbb_max_threads = 1;
#ifdef STAN_THREADS
  if (n_threads == 0) {
    tbb_max_threads = internal::get_num_threads();
  } else if (n_threads > 0) {
    tbb_max_threads = n_threads;
  } else if (n_threads == -1) {
    tbb_max_threads = std::thread::hardware_concurrency();
  } else {
    invalid_argument("init_threadpool_tbb(int)", "n_threads", n_threads,
                     "The number of threads is '",
                     "' but it must be positive or -1");
  }
#endif
  static tbb::global_control tbb_gc(
      tbb::global_control::max_allowed_parallelism, tbb_max_threads);
  static tbb::task_arena tbb_arena(tbb_max_threads, 1);
  tbb_arena.initialize();

  return tbb_arena;
}
#else
/**
 * Initialize the Intel TBB threadpool and global scheduler through
 * the tbb::task_scheduler_init object. In case an instance of the
 * tbb::task_scheduler_object has been instantiated prior to calling
 * this function, then any subsequent initialization is ignored by the
 * Intel TBB.
 *
 * The maximal number of threads is read from the environment variable
 * STAN_NUM_THREADS using internal::get_num_threads. See conventions
 * of get_num_threads. The TBB scheduler will be activated by calling
 * this function.
 *
 * The function returns a reference to the static
 * tbb::task_scheduler_init instance.
 *
 * @param n_threads The maximum number of threads available to the tbb. If not
 *  set will search for the environment variable `STAN_NUM_THREADS`. A value of
 *  will assume all detectable threads are available for the process.
 * @return reference to the static tbb::task_scheduler_init
 * @throws std::runtime_error if n_threads (defaults to zero) is not provided
 * and the value of STAN_NUM_THREADS environment variable is invalid.
 */
inline tbb::task_scheduler_init& init_threadpool_tbb(int n_threads = 0) {
  int tbb_max_threads = 1;
#ifdef STAN_THREADS
  if (n_threads == 0) {
    tbb_max_threads = internal::get_num_threads();
  } else if (n_threads > 0) {
    tbb_max_threads = n_threads;
  } else if (n_threads == -1) {
    tbb_max_threads = std::thread::hardware_concurrency();
  } else {
    invalid_argument("init_threadpool_tbb(int)", "n_threads", n_threads,
                     "The number of threads is '",
                     "' but it must be positive or -1");
  }
#endif
  static tbb::task_scheduler_init tbb_scheduler(tbb_max_threads, 0);
  return tbb_scheduler;
}
#endif

}  // namespace math
}  // namespace stan

#endif
