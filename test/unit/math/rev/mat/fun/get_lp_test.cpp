#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/util.hpp>


TEST(mathMatrix,getLp) {
  using stan::math::accumulator;
  using stan::math::get_lp;
  using stan::math::var;

  var lp = 12.5;
  accumulator<var> lp_accum;
  EXPECT_FLOAT_EQ(12.5, get_lp(lp,lp_accum).val());

  lp_accum.add(2);
  lp_accum.add(3);
  EXPECT_FLOAT_EQ(17.5, get_lp(lp,lp_accum).val());
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::accumulator;
  using stan::math::get_lp;
  using stan::math::var;

  var lp = 12.5;
  accumulator<var> lp_accum;
  lp_accum.add(2);
  lp_accum.add(3);
  test::check_varis_on_stack(stan::math::get_lp(lp,lp_accum));
}
