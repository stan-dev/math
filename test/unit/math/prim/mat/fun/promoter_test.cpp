#include <stan/math/prim/mat.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/is_same.hpp>

TEST(MathFunctionsMatPromoter,int2double) {
  std::vector<int> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  stan::math::promoter<std::vector<int>,std::vector<double> > p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<std::vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsMatPromoter,double2double) {
  std::vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  stan::math::promoter<std::vector<double>,std::vector<double> > p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<std::vector<double>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsMatPromoter,double2var) {
  std::vector<double> from;
  from.push_back(1);
  from.push_back(2);
  from.push_back(3);
  stan::math::promoter<std::vector<double>,std::vector<stan::math::var> > p;
  typedef BOOST_TYPEOF(p.promote_to(from)) result_t;
  bool same = boost::is_same<std::vector<stan::math::var>, result_t>::value;
  EXPECT_TRUE(same);
}

