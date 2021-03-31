#include <stan/math/prim/core.hpp>
#include <gtest/gtest.h>


TEST(intel_tbb_new_init, set_threads) {
  int num_threads = std::thread::hardware_concurrency() - 1;
  auto& tbb_init = stan::math::init_threadpool_tbb(num_threads);
  EXPECT_TRUE(tbb_init.is_active());

  EXPECT_EQ(num_threads, tbb_init.max_concurrency());
}
