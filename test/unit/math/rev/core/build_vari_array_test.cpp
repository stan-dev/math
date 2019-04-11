#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <sstream>

TEST(AgradRev, build_vari_array) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using stan::math::vari;

  Matrix<var, Dynamic, Dynamic> mdd(3, 2);
  Matrix<var, 1, Dynamic> mvd(5);
  Matrix<var, Dynamic, 1> mdv(3);
  Matrix<var, 2, 3> mvv;

  mdd << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  mvd << 1.0, 2.0, 3.0, 4.0, 5.0;
  mdv << 1.0, 2.0, 3.0;
  mvv << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

  vari **vdd = build_vari_array(mdd);
  vari **vvd = build_vari_array(mvd);
  vari **vdv = build_vari_array(mdv);
  vari **vvv = build_vari_array(mvv);

  EXPECT_TRUE(stan::math::ChainableStack::instance().memalloc_.in_stack(vdd));
  EXPECT_TRUE(stan::math::ChainableStack::instance().memalloc_.in_stack(vvd));
  EXPECT_TRUE(stan::math::ChainableStack::instance().memalloc_.in_stack(vdv));
  EXPECT_TRUE(stan::math::ChainableStack::instance().memalloc_.in_stack(vvv));

  for (int i = 0; i < mdd.size(); ++i) {
    EXPECT_EQ(mdd.data()[i].vi_, vdd[i]);
    EXPECT_TRUE(
        stan::math::ChainableStack::instance().memalloc_.in_stack(vdd[i]));
  }
  for (int i = 0; i < mvd.size(); ++i) {
    EXPECT_EQ(mvd.data()[i].vi_, vvd[i]);
    EXPECT_TRUE(
        stan::math::ChainableStack::instance().memalloc_.in_stack(vvd[i]));
  }
  for (int i = 0; i < mdv.size(); ++i) {
    EXPECT_EQ(mdv.data()[i].vi_, vdv[i]);
    EXPECT_TRUE(
        stan::math::ChainableStack::instance().memalloc_.in_stack(vdv[i]));
  }
  for (int i = 0; i < mvv.size(); ++i) {
    EXPECT_EQ(mvv.data()[i].vi_, vvv[i]);
    EXPECT_TRUE(
        stan::math::ChainableStack::instance().memalloc_.in_stack(vvv[i]));
  }
}
