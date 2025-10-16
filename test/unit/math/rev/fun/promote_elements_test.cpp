#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <boost/typeof/typeof.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

using Eigen::Matrix;
using stan::math::promote_elements;
using stan::math::var;
using std::vector;

TEST(MathFunctionsScalPromote_Elements, double2var) {
  double from;
  promote_elements<var, double> p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<var, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromote_Elements, doubleVec2varVec) {
  vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promote_elements<vector<var>, vector<double> > p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<vector<var>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsMatPromote_Elements, doubleMat2varMat) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  promote_elements<Matrix<var, 2, 3>, Matrix<double, 2, 3> > p;
  typedef BOOST_TYPEOF(p.promote(m1)) result_t;
  bool same = std::is_same<Matrix<var, 2, 3>, result_t>::value;
  EXPECT_TRUE(same);
}
