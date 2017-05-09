#include <gtest/gtest.h>
#include <stan/math/fwd/core/fvar.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>

TEST(foo, bar) {
  using stan::math::detail::broadcast_array;
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> fv(1.0, 2.1);
  broadcast_array<fvar<var> > ba(fv);
  EXPECT_EQ(1.0, ba[0].val_.vi_->val_);
  EXPECT_EQ(1.0, ba[2].val_.vi_->val_);
}
