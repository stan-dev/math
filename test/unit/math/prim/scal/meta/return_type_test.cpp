#include <stan/math/prim/scal.hpp>
#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::return_type;
using stan::math::var;

TEST(MetaTraits, ReturnTypeDouble) {
  test::expect_same_type<double, return_type<double>::type>();
}

TEST(MetaTraits, ReturnTypeFloat) {
  test::expect_same_type<double, return_type<float>::type>();
}

TEST(MetaTraits, ReturnTypeInt) {
  test::expect_same_type<double, return_type<int>::type>();
}

TEST(MetaTraits, ReturnTypeDoubleMany) {
  test::expect_same_type<double, return_type<double,int,double,float,float,
    double,float,int,double>::type>();
}

TEST(MetaTraits, ReturnTypeTwentyOneParams) {
  test::expect_same_type<var, return_type<double,var,double,int,double,float,
     float,float,var,int,double,int,double,float,float,
     double,var,double,int,double,float>::type>();
}

TEST(MetaTraits, ReturnTypeVarLast) {
  test::expect_same_type<var, return_type<double,double,var>::type>();
}
