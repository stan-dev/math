#include <stan/math/prim/core.hpp>
#include <gtest/gtest.h>

TEST(intel_tbb_new_init, set_threads) {
  int num_threads = std::thread::hardware_concurrency();
  auto& tbb_init = stan::math::init_threadpool_tbb(num_threads - 1);
  EXPECT_TRUE(tbb_init.is_active());

#ifdef STAN_THREADS
  EXPECT_EQ(num_threads - 1, tbb_init.max_concurrency());
#else
  EXPECT_EQ(1, tbb_init.max_concurrency());
#endif
}
