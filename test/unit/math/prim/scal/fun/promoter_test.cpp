#include <stan/math/prim/scal.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/is_same.hpp>

using stan::math::promoter;
using stan::math::var;

TEST(MathFunctionsScalPromoter,int2double) {
  int from;
  promoter<int,double> p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsScalPromoter,double2double) {
  double from;
  promoter<double,double> p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsScalPromoter,double2var) {
  double from;
  promoter<double,var> p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<var, result_t>::value;
  EXPECT_TRUE(same);
}

