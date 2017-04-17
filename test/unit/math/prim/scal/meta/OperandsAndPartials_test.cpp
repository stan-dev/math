#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartials) {
  using stan::math::OperandsAndPartials;

  OperandsAndPartials<double> o1;
  OperandsAndPartials<double,double,double,double> o2;

  SUCCEED() << "Construction should be ok.";
}

// TEST(AgradPartialsVari, incr_deriv_double) {
//   using stan::VectorView;
//   using stan::math::incr_deriv;
//   using stan::is_vector;
//   using stan::is_constant_struct;

//   double a = 1;
  
//   VectorView<const double,stan::is_vector<double>::value,
//              stan::is_constant_struct<double>::value> d_a(stan::length(a) * 0);


  
//   double result = incr_deriv<VectorView<const double,
//                                            is_vector<double>::value,
//                                            is_constant_struct<double>::value>,
//                              double,double>().incr(d_a,a);
//   EXPECT_FLOAT_EQ(0, result);
// }


TEST(AgradPartialsVari, OperandsAndPartials_check_throw) {
  using stan::math::OperandsAndPartials;
  
  double d;
  
  OperandsAndPartials<> o1(d,d,d,d,d,d);
  EXPECT_THROW(o1.d_x1[0], std::logic_error);
  EXPECT_THROW(o1.d_x2[0], std::logic_error);
  EXPECT_THROW(o1.d_x3[0], std::logic_error);
  EXPECT_THROW(o1.d_x4[0], std::logic_error);
  EXPECT_THROW(o1.d_x5[0], std::logic_error);
  EXPECT_THROW(o1.d_x6[0], std::logic_error);
}
