#ifndef STAN_MATH_PRIM_SCAL_INIT_THREADPOOL_TBB_HPP
#define STAN_MATH_PRIM_SCAL_INIT_THREADPOOL_TBB_HPP

#ifdef STAN_THREADS

#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <boost/lexical_cast.hpp>

#include <tbb/task_scheduler_init.h>

#include <cstdlib>

namespace stan {
namespace math {

inline bool init_threadpool_tbb(bool use_env = true,
                                int max_threads
                                = tbb::task_scheduler_init::automatic,
                                tbb::stack_size_type stack_size = 0) {
  int tbb_max_threads = max_threads;
  if (use_env) {
    const char* env_stan_num_threads = std::getenv("STAN_NUM_THREADS");
    if (env_stan_num_threads != nullptr) {
      try {
        const int env_num_threads
            = boost::lexical_cast<int>(env_stan_num_threads);
        if (env_num_threads > 0)
          tbb_max_threads = env_num_threads;
        else if (env_num_threads == -1)
          tbb_max_threads = std::thread::hardware_concurrency();
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
  }

  static tbb::task_scheduler_init tbb_scheduler(tbb_max_threads, stack_size);

  return tbb_scheduler.is_active();
}

}  // namespace math
}  // namespace stan

#endif

#endif
