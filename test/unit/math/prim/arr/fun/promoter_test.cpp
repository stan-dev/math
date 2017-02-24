#include <stan/math/prim/arr.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/is_same.hpp>

using std::vector;
using stan::math::promoter;
using stan::math::var;

TEST(MathFunctionsArrPromoter,intVec2doubleVec) {
  vector<int> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promoter<vector<int>,vector<double> > p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromoter,doubleVec2doubleVec) {
  vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promoter<vector<double>,vector<double> > p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromoter,doubleVec2varVec) {
  vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  promoter<vector<double>,vector<var> > p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<vector<var>, result_t>::value;
  EXPECT_TRUE(same);
}

