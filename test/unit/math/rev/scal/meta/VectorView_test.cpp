#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorView_var) {
  using stan::VectorView;
  using stan::math::var;
  
  var d(10);
  VectorView<var> dv(d);
  EXPECT_FLOAT_EQ(d.val(), dv[0].val());
  dv[1] = 2.0;
  EXPECT_FLOAT_EQ(2.0, dv[0].val());
  EXPECT_FLOAT_EQ(2.0, d.val());

  const var c(10);
  VectorView<const var> cv(c);
  EXPECT_FLOAT_EQ(c.val(), cv[0].val());
}
