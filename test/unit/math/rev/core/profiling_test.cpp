#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(Profiling, profile_double_basic) {
  using stan::math::profile;
  using stan::math::var;
  stan::math::profile_map profiles;
  double a = 3.0, b = 2.0, c;
  {
    profile<double> p1 = profile<double>("p1", profiles);
    c = a + b;
  }
  {
    profile<int> p1 = profile<int>("p1", profiles);
    c = a + b;
  }
  stan::math::profile_key key = {"p1", std::this_thread::get_id()};
  EXPECT_EQ(profiles[key].get_chain_stack_used(), 0);
  EXPECT_EQ(profiles[key].get_nochain_stack_used(), 0);
  EXPECT_FLOAT_EQ(profiles[key].get_rev_time(), 0.0);
  EXPECT_EQ(profiles[key].get_num_rev_passes(), 0);
  EXPECT_EQ(profiles[key].get_num_fwd_passes(), 2);
  EXPECT_TRUE(profiles[key].get_fwd_time() > 0.0);
}

TEST(Profiling, profiling_var_basic) {
  using stan::math::profile;
  using stan::math::var;
  stan::math::profile_map profiles;
  var c;
  {
    profile<var> t1 = profile<var>("t1", profiles);
    var a = 2.0;
    var b = 3.0;
    c = a + b;
  }
  c.grad();
  stan::math::recover_memory();
  stan::math::profile_key key = {"t1", std::this_thread::get_id()};
  EXPECT_EQ(profiles[key].get_chain_stack_used(), 1);
  EXPECT_EQ(profiles[key].get_nochain_stack_used(), 2);
  EXPECT_EQ(profiles[key].get_num_rev_passes(), 1);
  EXPECT_EQ(profiles[key].get_num_fwd_passes(), 1);
  EXPECT_TRUE(profiles[key].get_fwd_time() > 0.0); 
  EXPECT_TRUE(profiles[key].get_rev_time() > 0.0);
}

