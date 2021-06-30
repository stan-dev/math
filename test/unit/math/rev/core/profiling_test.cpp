#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(Profiling, double_basic) {
  using stan::math::profile;
  using stan::math::var;
  stan::math::profile_map profiles;
  double a = 3.0, b = 2.0, c;
  {
    profile<double> p1("p1", profiles);
    c = log(exp(a)) * log(exp(b));
  }
  {
    profile<int> p1("p1", profiles);
    c = log(exp(a)) * log(exp(b));
  }

  stan::math::profile_key key = {"p1", std::this_thread::get_id()};
  EXPECT_NEAR(c, 6.0, 1E-8);
  EXPECT_EQ(profiles[key].get_chain_stack_used(), 0);
  EXPECT_EQ(profiles[key].get_nochain_stack_used(), 0);
  EXPECT_FLOAT_EQ(profiles[key].get_rev_time(), 0.0);
  EXPECT_EQ(profiles[key].get_num_rev_passes(), 0);
  EXPECT_EQ(profiles[key].get_num_fwd_passes(), 2);
  EXPECT_EQ(profiles[key].get_num_no_AD_fwd_passes(), 2);
  EXPECT_EQ(profiles[key].get_num_AD_fwd_passes(), 0);
  EXPECT_TRUE(profiles[key].get_fwd_time() > 0.0);
}

TEST(Profiling, var_basic) {
  using stan::math::profile;
  using stan::math::var;
  stan::math::profile_map profiles;
  var c;
  {
    profile<var> t1("t1", profiles);
    var a = 2.0;
    var b = 3.0;
    c = a + b;
  }
  c.grad();
  stan::math::recover_memory();
  stan::math::profile_key key = {"t1", std::this_thread::get_id()};
  EXPECT_EQ(profiles[key].get_chain_stack_used(), 1);
  EXPECT_EQ(profiles[key].get_nochain_stack_used(), 2);
  EXPECT_EQ(profiles[key].get_num_no_AD_fwd_passes(), 0);
  EXPECT_EQ(profiles[key].get_num_AD_fwd_passes(), 1);
  EXPECT_EQ(profiles[key].get_num_rev_passes(), 1);
  EXPECT_EQ(profiles[key].get_num_fwd_passes(), 1);
  EXPECT_TRUE(profiles[key].get_fwd_time() > 0.0);
  EXPECT_TRUE(profiles[key].get_rev_time() > 0.0);
}

TEST(Profiling, var_exception) {
  using stan::math::profile;
  using stan::math::var;
  stan::math::profile_map profiles;
  var c;
  try {
    profile<var> t1("t1", profiles);
    var a = 2.0;
    var b = 3.0;
    c = a + b;
    throw std::domain_error("error");
    c.grad();
  } catch (const std::exception& e) {
  }
  stan::math::recover_memory();
  stan::math::profile_key key_t1 = {"t1", std::this_thread::get_id()};
  EXPECT_EQ(profiles[key_t1].get_chain_stack_used(), 1);
  EXPECT_EQ(profiles[key_t1].get_nochain_stack_used(), 2);
  EXPECT_EQ(profiles[key_t1].get_num_rev_passes(), 0);
  EXPECT_EQ(profiles[key_t1].get_num_fwd_passes(), 1);
  EXPECT_EQ(profiles[key_t1].get_num_no_AD_fwd_passes(), 0);
  EXPECT_EQ(profiles[key_t1].get_num_AD_fwd_passes(), 1);
  EXPECT_EQ(profiles[key_t1].get_num_rev_passes(), 0);
  EXPECT_EQ(profiles[key_t1].get_num_fwd_passes(), 1);
  EXPECT_TRUE(profiles[key_t1].get_fwd_time() > 0.0);
  EXPECT_TRUE(profiles[key_t1].get_rev_time() == 0.0);
}

TEST(Profiling, var_loop) {
  using stan::math::profile;
  using stan::math::var;
  stan::math::profile_map profiles;
  var c = 1.0;
  int N = 100;
  for (int i = 0; i < N; i++) {
    profile<var> t1("t1", profiles);
    var a = i;
    var b = 2.0;
    c = c + a;
  }
  c.grad();
  stan::math::recover_memory();
  stan::math::profile_key key_t1 = {"t1", std::this_thread::get_id()};
  EXPECT_EQ(profiles[key_t1].get_chain_stack_used(), N);
  EXPECT_EQ(profiles[key_t1].get_nochain_stack_used(), 2 * N);
  EXPECT_EQ(profiles[key_t1].get_num_rev_passes(), N);
  EXPECT_EQ(profiles[key_t1].get_num_fwd_passes(), N);
  EXPECT_EQ(profiles[key_t1].get_num_no_AD_fwd_passes(), 0);
  EXPECT_EQ(profiles[key_t1].get_num_AD_fwd_passes(), N);
  EXPECT_TRUE(profiles[key_t1].get_fwd_time() > 0.0);
  EXPECT_TRUE(profiles[key_t1].get_rev_time() > 0.0);
}

TEST(Profiling, duplicate_active_profile) {
  using stan::math::profile;
  using stan::math::var;
  stan::math::profile_map profiles;
  profile<var> t1("t1", profiles);
  EXPECT_THROW(profile<var>("t1", profiles), std::runtime_error);
}
