#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorView_vector_var) {
  using stan::VectorView;
  using std::vector;
  using stan::math::var;
  
  vector<var> x(10);
  for (size_t n = 0; n < 10; ++n) 
    x[n] = n;
  VectorView<vector<var> > xv(x);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(x[n].val(), xv[n].val());
  for (size_t n = 0; n < 10; ++n)
    xv[n] = 10+n;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(x[n].val(), xv[n].val());
    EXPECT_FLOAT_EQ(10+n, xv[n].val());
  }

  const vector<var> y(x);
  VectorView<const vector<var> > yv(y);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(y[n].val(), yv[n].val());
}
