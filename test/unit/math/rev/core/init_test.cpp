#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>

#include <thread>

TEST(init, thread_initialize) {
  // the main thread must be initialized by the time this code is
  // reached
  EXPECT_TRUE(stan::math::ChainableStack::instance_ != nullptr);

  stan::math::ChainableStack::AutodiffStackStorage& main_ad_stack
      = stan::math::ChainableStack::instance();

  EXPECT_TRUE(&main_ad_stack == stan::math::init());

#ifdef STAN_THREADS
  auto thread_tester = [&]() -> void {
    EXPECT_TRUE(stan::math::ChainableStack::instance_ == nullptr);
    stan::math::init();
    EXPECT_TRUE(stan::math::ChainableStack::instance_ != nullptr);
    EXPECT_TRUE(stan::math::ChainableStack::instance_ != &main_ad_stack);
  };
#else
  auto thread_tester = [&]() -> void {
    EXPECT_TRUE(stan::math::ChainableStack::instance_ != nullptr);
    stan::math::init();
    EXPECT_TRUE(stan::math::ChainableStack::instance_ != nullptr);
    EXPECT_TRUE(stan::math::ChainableStack::instance_ == &main_ad_stack);
  };
#endif
  std::thread other_work(thread_tester);

  other_work.join();
}

TEST(init, thread_instances) {
  // place a var on the stack such that a fresh stack in another
  // thread will be different at initialization (if STAN_THREADS is
  // set)
  stan::math::var a = 1;

  stan::math::ChainableStack::AutodiffStackStorage& main_ad_stack
      = stan::math::ChainableStack::instance();

#ifdef STAN_THREADS
  auto thread_tester = [&]() -> void {
    stan::math::init();
    EXPECT_TRUE(main_ad_stack.var_stack_.size()
                > stan::math::ChainableStack::instance().var_stack_.size());
  };
#else
  auto thread_tester = [&]() -> void {
    stan::math::init();
    EXPECT_TRUE(main_ad_stack.var_stack_.size()
                == stan::math::ChainableStack::instance().var_stack_.size());
  };
#endif

  std::thread other_work(thread_tester);

  other_work.join();
}
