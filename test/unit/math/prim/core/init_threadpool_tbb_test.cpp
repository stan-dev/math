#include <stan/math/prim/core.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>

TEST(intel_tbb_new_init, check_status) {
  set_n_threads(-1);
  auto& tbb_init = stan::math::init_threadpool_tbb();
  EXPECT_TRUE(tbb_init.is_active());

  EXPECT_EQ(std::thread::hardware_concurrency(), tbb_init.max_concurrency());

  auto& tbb_reinit = stan::math::init_threadpool_tbb();
  EXPECT_TRUE(tbb_init.is_active());

  tbb_init.terminate();
  EXPECT_FALSE(tbb_init.is_active());
}
