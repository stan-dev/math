#include <stan/math/prim/core.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>

TEST(intel_tbb_init, set_max_threads) {
  tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb(10);
  EXPECT_TRUE(tbb_init.is_active());

  EXPECT_EQ(10, tbb::this_task_arena::max_concurrency());

  tbb_init.terminate();
  EXPECT_FALSE(tbb_init.is_active());

  tbb_init = stan::math::init_threadpool_tbb(4);
  EXPECT_TRUE(tbb_init.is_active());

  EXPECT_EQ(4, tbb::this_task_arena::max_concurrency());

  // tbb_init_2.terminate();
  // EXPECT_FALSE(tbb_init_2.is_active());

  // tbb::task_scheduler_init& tbb_init_3 = stan::math::init_threadpool_tbb(1);
  // EXPECT_TRUE(tbb_init_3.is_active());

  // EXPECT_EQ(1, tbb::this_task_arena::max_concurrency());

  // tbb_init_3.terminate();
  // EXPECT_FALSE(tbb_init_3.is_active());
}

TEST(intel_tbb_init, check_status) {
  tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb(-1);
  EXPECT_TRUE(tbb_init.is_active());

  EXPECT_EQ(std::thread::hardware_concurrency(),
            tbb::this_task_arena::max_concurrency());

  tbb::task_scheduler_init& tbb_reinit = stan::math::init_threadpool_tbb(-1);
  EXPECT_TRUE(tbb_init.is_active());

  tbb_init.terminate();
  EXPECT_FALSE(tbb_init.is_active());
}

// TEST(intel_tbb_init, incorrect_values) {
//   EXPECT_THROW_MSG(tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb(0),
//                    std::invalid_argument, "positive number or -1");

//   EXPECT_THROW_MSG(tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb(-500),
//                    std::invalid_argument, "positive number or -1");

//   EXPECT_THROW_MSG(tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb(-3),
//                    std::invalid_argument, "must be positive or -1");
// }

// // tests with deprecated STAN_NUM_THREADS env. variable
// TEST(intel_tbb_init, set_max_threads_deprecated) {
//   set_n_threads("10");
//   tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb(-1);
//   EXPECT_TRUE(tbb_init.is_active());
//   EXPECT_EQ(10, tbb::this_task_arena::max_concurrency());

//   tbb_init.terminate();
//   EXPECT_FALSE(tbb_init.is_active());

//   set_n_threads("4");
//   tbb::task_scheduler_init& tbb_init_2 = stan::math::init_threadpool_tbb(-1);
//   EXPECT_TRUE(tbb_init_2.is_active());
//   EXPECT_EQ(4, tbb::this_task_arena::max_concurrency());

//   tbb_init_2.terminate();
//   EXPECT_FALSE(tbb_init_2.is_active());

//   set_n_threads("-1");
//   tbb::task_scheduler_init& tbb_init_3 = stan::math::init_threadpool_tbb(-1);
//   EXPECT_TRUE(tbb_init_3.is_active());
//   EXPECT_EQ(std::thread::hardware_concurrency(), tbb::this_task_arena::max_concurrency());

//   tbb_init_3.terminate();
//   EXPECT_FALSE(tbb_init_3.is_active());
// }

// TEST(intel_tbb_init, incorrect_values_deprecated) {
//   set_n_threads("abc");
//   EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(-1),
//                    std::invalid_argument, "positive number or -1");

//   set_n_threads("1c");
//   EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(-1),
//                    std::invalid_argument, "positive number or -1");

//   set_n_threads("-2");
//   EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(-1),
//                    std::invalid_argument, "positive number or -1");

//   set_n_threads("0");
//   EXPECT_THROW_MSG(stan::math::init_threadpool_tbb(-1),
//                    std::invalid_argument, "positive number or -1");
// }


// TEST(intel_tbb_init, check_status_deprecated) {
//   set_n_threads(-1);
//   tbb::task_scheduler_init& tbb_init = stan::math::init_threadpool_tbb();
//   EXPECT_TRUE(tbb_init.is_active());

//   EXPECT_EQ(std::thread::hardware_concurrency(),
//             tbb::this_task_arena::max_concurrency());

//   tbb::task_scheduler_init& tbb_reinit = stan::math::init_threadpool_tbb();
//   EXPECT_TRUE(tbb_init.is_active());

//   tbb_init.terminate();
//   EXPECT_FALSE(tbb_init.is_active());
// }
