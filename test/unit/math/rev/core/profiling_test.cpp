#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(Prolfiling, chain_stack_test) {
  using stan::math::var;
  using stan::math::profile;
  stan::math::profile_map profiles;
  var a, b, c;
  {
        profile<var> t1 = profile<var>("t1", profiles);
        a = 2.0;
        b = 3.0;
        c = a + b;
  }
  c.grad();
  stan::math::profile_key key = {"t1", std::this_thread::get_id()};
  EXPECT_EQ(profiles[key].chain_stack_size_sum, 1);
  EXPECT_EQ(profiles[key].nochain_stack_size_sum, 1);
  stan::math::recover_memory();
}
