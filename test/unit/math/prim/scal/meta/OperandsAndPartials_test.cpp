#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, incr_deriv_double) {
  using stan::VectorView;
  using stan::math::incr_deriv;
  using stan::is_vector;
  using stan::is_constant_struct;

  double a = 1;
  
  VectorView<const double,stan::is_vector<double>::value,
             stan::is_constant_struct<double>::value> d_a(stan::length(a) * 0);


  
  double result = incr_deriv<VectorView<const double,
                                           is_vector<double>::value,
                                           is_constant_struct<double>::value>,
                             double,double>().incr(d_a,a);
  EXPECT_FLOAT_EQ(0, result);
}
