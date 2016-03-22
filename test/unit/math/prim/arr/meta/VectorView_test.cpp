#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorView_vector_double) {
  using stan::VectorView;
  using std::vector;
  
  vector<double> x(10);
  for (size_t n = 0; n < 10; ++n) 
    x[n] = n;
  VectorView<vector<double> > xv(x);
  for (size_t n = 0; n < 10; ++n)
    EXPECT_FLOAT_EQ(x[n], xv[n]);
  for (size_t n = 0; n < 10; ++n)
    xv[n] = 10+n;
  for (size_t n = 0; n < 10; ++n) {
    EXPECT_FLOAT_EQ(x[n], xv[n]);
    EXPECT_FLOAT_EQ(10+n, xv[n]);
  }

  const vector<double> y(x);
  VectorView<const vector<double> > yv(y);
  // for (size_t n = 0; n < 10; ++n)
  //   EXPECT_FLOAT_EQ(y[n], yv[n]);
}

TEST(MetaTraits,VectorView) {
  using stan::VectorView;
  using std::vector;

  vector<double> sv;
  sv.push_back(1.0);
  sv.push_back(4.0);
  sv.push_back(9.0);
  // VectorView<const vector<double> > sv_VectorView(sv);
  // EXPECT_FLOAT_EQ(1.0,sv_VectorView[0]);
  // EXPECT_FLOAT_EQ(4.0,sv_VectorView[1]);
  // EXPECT_FLOAT_EQ(9.0,sv_VectorView[2]);
}
