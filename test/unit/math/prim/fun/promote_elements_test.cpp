#include <stan/math/prim.hpp>
#include <boost/typeof/typeof.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

TEST(MathFunctionsScalPromote_Elements, int2double) {
  using Eigen::Matrix;
  using stan::math::promote_elements;
  using std::vector;
  int from;
  promote_elements<double, int> p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsScalPromote_Elements, double2double) {
  using Eigen::Matrix;
  using stan::math::promote_elements;
  using std::vector;
  double from;
  promote_elements<double, double> p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromote_Elements, intVec2doubleVec) {
  using Eigen::Matrix;
  using stan::math::promote_elements;
  using std::vector;
  vector<int> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promote_elements<vector<double>, vector<int> > p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromote_Elements, doubleVec2doubleVec) {
  using Eigen::Matrix;
  using stan::math::promote_elements;
  using std::vector;
  vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promote_elements<vector<double>, vector<double> > p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsMatPromote_Elements, doubleMat2doubleMat) {
  using Eigen::Matrix;
  using stan::math::promote_elements;
  using std::vector;
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  promote_elements<Matrix<double, 2, 3>, Matrix<double, 2, 3> > p;
  typedef BOOST_TYPEOF(p.promote(m1)) result_t;
  bool same = std::is_same<Matrix<double, 2, 3>, result_t>::value;
  EXPECT_TRUE(same);
}
