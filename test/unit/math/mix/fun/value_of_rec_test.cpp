#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>

TEST(AgradMix, value_of_rec) {
  using stan::math::fvar;
  using stan::math::value_of_rec;
  using stan::math::var;

  fvar<var> fv_a(5.0);
  fvar<fvar<var> > ffv_a(5.0);
  fvar<fvar<fvar<fvar<fvar<var> > > > > fffffv_a(5.0);

  EXPECT_FLOAT_EQ(5.0, value_of_rec(fv_a));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(ffv_a));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(fffffv_a));
}
#include <stan/math/mix/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixMixArr, value_of_rec) {
  using stan::math::fvar;
  using stan::math::value_of_rec;
  using stan::math::var;
  using std::vector;

  vector<fvar<fvar<var> > > a;
  for (size_t i = 0; i < 10; ++i)
    a.push_back(fvar<fvar<var> >(i + 1));

  vector<fvar<fvar<var> > > b;
  for (size_t i = 10; i < 15; ++i)
    b.push_back(fvar<fvar<var> >(i + 1));

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val_.val_.val(), d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val_.val_.val(), d_a[i]);
}
