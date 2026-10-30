#include <stan/math/prim/core.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(intel_tbb_new_init, set_threads) {
  using control = tbb::global_control;
  const auto previous = control::active_value(control::max_allowed_parallelism);
  for (int num_threads : {1, 2, 4}) {
    {
      auto state = stan::math::internal::make_threadpool_tbb(num_threads);
      EXPECT_TRUE(state.arena.is_active());
      EXPECT_EQ(num_threads, state.arena.max_concurrency());
      EXPECT_EQ(std::min(previous, static_cast<std::size_t>(num_threads)),
                control::active_value(control::max_allowed_parallelism));
      state.arena.execute([num_threads]() {
        EXPECT_EQ(num_threads, tbb::this_task_arena::max_concurrency());
      });
    }
    EXPECT_EQ(previous,
              control::active_value(control::max_allowed_parallelism));
  }
}
