#include <stan/math/prim/core.hpp>

#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>
#include <tbb/task_arena.h>
#include <gtest/gtest.h>

#ifdef TBB_INTERFACE_NEW

#include <tbb/global_control.h>

TEST(intel_tbb_new_late_init, check_status) {
  const int num_threads = tbb::global_control::max_allowed_parallelism;

  if (num_threads > 1) {
    set_n_threads(num_threads - 1);
    tbb::task_arena& tbb_arena = stan::math::init_threadpool_tbb();
    tbb_arena.execute([&]() {
      EXPECT_TRUE(tbb_arena.is_active());

      // STAN_NUM_THREADS is not being honored if we have first
      // initialized the TBB scheduler outside of init_threadpool_tbb
#ifdef STAN_THREADS
      EXPECT_EQ(num_threads, tbb::this_task_arena::max_concurrency());
#else
      EXPECT_EQ(num_threads, 1);
#endif
    });
  }
}

#else

#include <tbb/task_scheduler_init.h>

TEST(intel_tbb_late_init, check_status) {
  const int num_threads = tbb::task_scheduler_init::default_num_threads();
  tbb::task_scheduler_init tbb_scheduler;

  if (num_threads > 1) {
    set_n_threads(num_threads - 1);
    tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb();
    EXPECT_TRUE(tbb_init.is_active());

    // STAN_NUM_THREADS is not being honored if we have first
    // initialized the TBB scheduler outside of init_threadpool_tbb
#ifdef STAN_THREADS
    EXPECT_EQ(num_threads, tbb::this_task_arena::max_concurrency());
#else
    EXPECT_EQ(num_threads, 1);
#endif
  }
}

#endif
