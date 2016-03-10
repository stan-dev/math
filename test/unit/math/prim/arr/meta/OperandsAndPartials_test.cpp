#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

// TEST(AgradPartialsVari, incr_deriv_vec_double) {
//   using stan::VectorView;
//   using stan::math::incr_deriv;
//   using stan::is_vector;
//   using stan::is_constant_struct;

//   std::vector<double> b;
//   b.push_back(1);
//   b.push_back(1);
  
//   VectorView<const std::vector<double>,
//              stan::is_vector<std::vector<double> >::value,
//              stan::is_constant_struct<std::vector<double> >::value> 
//     d_b(stan::length(b) * 0);

//   double result = incr_deriv<VectorView<const std::vector<double>,
//                                         is_vector<std::vector<double> >::value,
//                                         is_constant_struct<std::vector<double> >::value>,
//                              std::vector<double>,double>().incr(d_b,b);
  
//   EXPECT_FLOAT_EQ(0, result);
// }

TEST(AgradPartialsVari, OperandsAndPartials_check_throw) {
  using stan::math::OperandsAndPartials;
  using std::vector;
  
  vector<double> D;
  
  OperandsAndPartials<vector<double>,vector<double>,vector<double>,
                      vector<double>,vector<double>,vector<double> > o3(D,D,D,D,D,D);
  EXPECT_THROW(o3.d_x1[0], std::logic_error);
  EXPECT_THROW(o3.d_x2[0], std::logic_error);
  EXPECT_THROW(o3.d_x3[0], std::logic_error);
  EXPECT_THROW(o3.d_x4[0], std::logic_error);
  EXPECT_THROW(o3.d_x5[0], std::logic_error);
  EXPECT_THROW(o3.d_x6[0], std::logic_error);
}
