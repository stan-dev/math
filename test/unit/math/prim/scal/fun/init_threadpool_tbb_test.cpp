#include <stan/math/prim/core.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/mat/functor/utils_threads.hpp>

#ifdef STAN_THREADS
#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>
#endif

TEST(intel_tbb_init, check_status) {
  set_n_threads(-1);
  const bool tbb_init = stan::math::init_threadpool_tbb();
#ifdef STAN_THREADS
  EXPECT_TRUE(tbb_init);

  EXPECT_EQ(tbb::task_scheduler_init::default_num_threads(),
            tbb::this_task_arena::max_concurrency());
#else
  EXPECT_FALSE(tbb_init);
#endif
}

// TODO(wds15):
// - test stack_size argument being set
// - test that active status is false if we init before another
// task_scheduler_init
