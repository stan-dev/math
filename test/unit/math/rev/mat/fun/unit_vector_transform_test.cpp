#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/jacobian.hpp>
#include <test/unit/math/rev/mat/util.hpp>

using Eigen::Matrix;
using Eigen::Dynamic;

TEST(probTransform,unit_vector_jacobian) {
  using stan::math::var;
  using std::vector;
  var a = 2.0;
  var b = 3.0;
  var c = -1.0;
  
  Matrix<var,Dynamic,1> y(3);
  y << a, b, c;
  
  var lp(0);
  Matrix<var,Dynamic,1> x 
    = stan::math::unit_vector_constrain(y,lp);
  const var r2 = stan::math::dot_self(y);

  vector<var> indeps;
  indeps.push_back(a);
  indeps.push_back(b);
  indeps.push_back(c);

  vector<var> deps;
  deps.push_back(x(0));
  deps.push_back(x(1));
  deps.push_back(r2);
  
  vector<vector<double> > jacobian;
  stan::math::jacobian(deps,indeps,jacobian);

  Matrix<double,Dynamic,Dynamic> J(3,3);
  for (int m = 0; m < 3; ++m) {
    for (int n = 0; n < 3; ++n) {
      J(m,n) = jacobian[m][n];
    }
  }
  
  double det_J = J.determinant();

  EXPECT_FLOAT_EQ(1.0 / det_J, lp.val())
    << "J = " << std::endl
    << J << std::endl
    << "det_J = " << det_J << std::endl
    << "x = " << x.transpose();
}

TEST(probTransform, check_varis_on_stack) {
  using stan::math::var;
  Matrix<var,Dynamic,1> y(3);
  y << 2, 3, -1;
  
  var lp(0);
  test::check_varis_on_stack(stan::math::unit_vector_constrain(y, lp));
  test::check_varis_on_stack(stan::math::unit_vector_constrain(y));
}
