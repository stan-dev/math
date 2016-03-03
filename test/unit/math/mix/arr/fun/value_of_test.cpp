#include <stan/math/mix/arr.hpp>
#include <gtest/gtest.h>

template <typename T>
void fill(const std::vector<double>& contents,
          std::vector<T>& x){
  x.assign(contents.size(), T());
  for (size_t i = 0; i < contents.size(); ++i)
    x[i] = T(contents[i]);
}

TEST(MathMatrix,value_of) {
  using stan::math::value_of;
  using std::vector;
  using stan::math::fvar;
  using stan::math::var;

  vector<double> a_vals;

  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;

  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  vector<fvar<fvar<var> > > a;
  fill(a_vals, a);
  vector<fvar<fvar<var> > > b;
  fill(b_vals, b);

  vector<fvar<var> > d_a = value_of(a);
  vector<fvar<var> > d_b = value_of(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val_.val_.val(), d_b[i].val_.val());

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val_.val_.val(), d_a[i].val_.val());

}
