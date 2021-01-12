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
  EXPECT_EQ(profiles[key].chain_stack_size_sum, 0);
  EXPECT_EQ(profiles[key].nochain_stack_size_sum, 0);
  EXPECT_FLOAT_EQ(profiles[key].rev_pass_time, 0.0);
  EXPECT_EQ(profiles[key].n_rev_pass, 0);
  EXPECT_EQ(profiles[key].n_fwd_pass, 2);
  EXPECT_TRUE(profiles[key].fwd_pass_time > 0.0);
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
  EXPECT_EQ(profiles[key].chain_stack_size_sum, 1);
  EXPECT_EQ(profiles[key].nochain_stack_size_sum, 2);
  EXPECT_EQ(profiles[key].n_rev_pass, 1);
  EXPECT_EQ(profiles[key].n_fwd_pass, 1);
  EXPECT_TRUE(profiles[key].fwd_pass_time > 0.0);
  EXPECT_TRUE(profiles[key].rev_pass_time > 0.0);  
}

