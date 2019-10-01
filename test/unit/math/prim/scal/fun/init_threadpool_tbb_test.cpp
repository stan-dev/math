#include <stan/math/prim/scal.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/mat/functor/utils_threads.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>

TEST(intel_tbb_init, check_status) {
  set_n_threads(-1);
  const bool tbb_init = stan::math::init_threadpool_tbb();
#ifdef STAN_THREADS
  EXPECT_TRUE(tbb_init);

  EXPECT_EQ(std::thread::hardware_concurrency(),
            tbb::this_task_arena::max_concurrency());
#else
  EXPECT_FALSE(tbb_init);
#endif
}

// TODO:
// - test stack_size argument being set
// - test that active status is false if we init before another
// task_scheduler_init
