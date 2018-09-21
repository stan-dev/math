#include <stan/math/mix/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrix, value_of_rec) {
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
