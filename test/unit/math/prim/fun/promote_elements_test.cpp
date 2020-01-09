#include <stan/math/prim/scal.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <type_traits>

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
