#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(AgradPartialsVari, OperandsAndPartials) {
  using stan::math::OperandsAndPartials;
  using stan::math::var;

  std::vector<double> d_vec(4);
  OperandsAndPartials<std::vector<double> > o3(d_vec);
  std::vector<var> v_vec;
  v_vec.push_back(var(0.0));
  v_vec.push_back(var(1.0));
  v_vec.push_back(var(2.0));
  v_vec.push_back(var(3.0));
  
  std::vector<double> grad;

  OperandsAndPartials<std::vector<var> > o4(v_vec);
  o4.d_x1[0] = 10.0;
  o4.d_x1[1] = 20.0;
  o4.d_x1[2] = 30.0;
  o4.d_x1[3] = 40.0;
  
  var v = o4.value(10.0);
  v.grad(v_vec, grad);
  EXPECT_EQ(4U, o4.nvaris);
  EXPECT_FLOAT_EQ(10.0, v.val());
  EXPECT_FLOAT_EQ(10.0, grad[0]);
  EXPECT_FLOAT_EQ(20.0, grad[1]);
  EXPECT_FLOAT_EQ(30.0, grad[2]);
  EXPECT_FLOAT_EQ(40.0, grad[3]);
}

TEST(AgradPartialsVari, OperandsAndPartials_check_throw) {
  using stan::math::OperandsAndPartials;
  using stan::math::var;
  using std::vector;
  
  vector<var> V;

  OperandsAndPartials<vector<var>,vector<var>,vector<var>,
                      vector<var>,vector<var>,vector<var> > o4(V,V,V,V,V,V);
  EXPECT_NO_THROW(o4.d_x1[0]);
  EXPECT_NO_THROW(o4.d_x2[0]);
  EXPECT_NO_THROW(o4.d_x3[0]);
  EXPECT_NO_THROW(o4.d_x4[0]);
  EXPECT_NO_THROW(o4.d_x5[0]);
  EXPECT_NO_THROW(o4.d_x6[0]);
}
