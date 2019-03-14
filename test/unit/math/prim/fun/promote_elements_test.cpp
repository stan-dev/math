
#include <stan/math/prim.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <type_traits>
#include <vector>
















using stan::math::promote_elements;
using stan::math::var;

TEST(MathFunctionsScalPromote_Elements, int2double) {
  int from;
  promote_elements<double, int> p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsScalPromote_Elements, double2double) {
  double from;
  promote_elements<double, double> p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsScalPromote_Elements, double2var) {
  double from;
  promote_elements<var, double> p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<var, result_t>::value;
  EXPECT_TRUE(same);
}







using stan::math::promote_elements;
using stan::math::var;
using std::vector;

TEST(MathFunctionsArrPromote_Elements_arr, intVec2doubleVec) {
  vector<int> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promote_elements<vector<double>, vector<int> > p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromote_Elements_arr, doubleVec2doubleVec) {
  vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promote_elements<vector<double>, vector<double> > p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromote_Elements_arr, doubleVec2varVec) {
  vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promote_elements<vector<var>, vector<double> > p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<vector<var>, result_t>::value;
  EXPECT_TRUE(same);
}






using Eigen::Matrix;
using stan::math::promote_elements;
using stan::math::var;

TEST(MathFunctionsMatPromote_Elements_mat, doubleMat2doubleMat) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  promote_elements<Matrix<double, 2, 3>, Matrix<double, 2, 3> > p;
  typedef BOOST_TYPEOF(p.promote(m1)) result_t;
  bool same = std::is_same<Matrix<double, 2, 3>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsMatPromote_Elements_mat, doubleMat2varMat) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  promote_elements<Matrix<var, 2, 3>, Matrix<double, 2, 3> > p;
  typedef BOOST_TYPEOF(p.promote(m1)) result_t;
  bool same = std::is_same<Matrix<var, 2, 3>, result_t>::value;
  EXPECT_TRUE(same);
}
