#include <stan/math/mix/arr.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix,value_of) {
  using stan::math::value_of;
  using std::vector;
  using stan::math::fvar;
  using stan::math::var;

  vector<fvar<fvar<var> > > a(10);
  for (size_t i = 0; i < 10; ++i)
    a[i] = fvar<fvar<var> >(i);
  vector<fvar<fvar<var> > > b(5);
  for (size_t i = 0; i < 5; ++i)
    b[i] = fvar<fvar<var> >(10 + i);

  vector<fvar<var> > d_a = value_of(a);
  vector<fvar<var> > d_b = value_of(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val_.val_.val(), d_b[i].val_.val());

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val_.val_.val(), d_a[i].val_.val());

}
