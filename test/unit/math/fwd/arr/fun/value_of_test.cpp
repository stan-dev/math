#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <gtest/gtest.h>

template <typename T>
void fill(const std::vector<double>& contents,
          std::vector<T>& x){
  size_t ij = 0;
  x.assign(contents.size(), T());
  for (int i = 0; i < contents.size(); ++i)
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

  {
    vector<fvar<double> > a;
    fill(a_vals, a);
    vector<fvar<double> > b;
    fill(b_vals, b);

    vector<double> d_a = value_of(a);
    vector<double> d_b = value_of(b);

    for (int i = 0; i < 5; ++i)
      EXPECT_FLOAT_EQ(b[i].val_, d_b[i]);

    for (int i = 0; i < 10; ++i)
      EXPECT_FLOAT_EQ(a[i].val_, d_a[i]);
  }

  {
    vector<fvar<fvar<double> > > a;
    fill(a_vals, a);
    vector<fvar<fvar<double> > > b;
    fill(b_vals, b);

    vector<fvar<double> > d_a = value_of(a);
    vector<fvar<double> > d_b = value_of(b);

    for (int i = 0; i < 5; ++i)
      EXPECT_FLOAT_EQ(b[i].val_.val_, d_b[i].val_);

    for (int i = 0; i < 10; ++i)
      EXPECT_FLOAT_EQ(a[i].val_.val_, d_a[i].val_);
  }

  {
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

}
