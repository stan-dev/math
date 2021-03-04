#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, set_zero_all_adjoints_nested_outside) {
  stan::math::var chaining = new stan::math::vari(1.0, true);
  chaining.adj() = 2.0;

  stan::math::start_nested();
  EXPECT_FLOAT_EQ(chaining.adj(), 2.0);
  stan::math::set_zero_all_adjoints_nested();
  EXPECT_FLOAT_EQ(chaining.adj(), 2.0);

  stan::math::recover_memory_nested();

  stan::math::var non_chaining = new stan::math::vari(1.0, false);
  non_chaining.adj() = 2.0;

  stan::math::start_nested();
  EXPECT_FLOAT_EQ(non_chaining.adj(), 2.0);
  stan::math::set_zero_all_adjoints_nested();
  EXPECT_FLOAT_EQ(non_chaining.adj(), 2.0);
}
