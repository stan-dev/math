#include <stan/math/prim/arr.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <type_traits>
#include <vector>

using stan::math::promote_elements;
using stan::math::var;
using std::vector;

TEST(MathFunctionsArrPromote_Elements, intVec2doubleVec) {
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
  vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promote_elements<vector<double>, vector<double> > p;
  typedef BOOST_TYPEOF(p.promote(from)) result_t;
  bool same = std::is_same<vector<double>, result_t>::value;
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
