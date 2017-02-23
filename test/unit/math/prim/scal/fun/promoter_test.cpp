#include <stan/math/prim/scal.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/is_same.hpp>

TEST(MathFunctionsArrPromoter,int2double) {
  int from;
  stan::math::promoter<int,double> p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromoter,double2double) {
  double from;
  stan::math::promoter<double,double> p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<double, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsArrPromoter,double2var) {
  double from;
  stan::math::promoter<double,stan::math::var> p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<stan::math::var, result_t>::value;
  EXPECT_TRUE(same);
}

