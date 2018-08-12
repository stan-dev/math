#include <stan/math/prim/scal.hpp>
#include <stan/math/rev/scal.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::var;
using stan::return_type;

TEST(MetaTraits, ReturnTypeDouble) {
  test::expect_same_type<double, return_type<double>::type>();
}

TEST(MetaTraits, ReturnTypeFloat) {
  test::expect_same_type<double, return_type<float>::type>();
}

TEST(MetaTraits, ReturnTypeInt) {
  test::expect_same_type<double, return_type<int>::type>();
}

TEST(MetaTraits, ReturnTypeVar) {
  test::expect_same_type<var, return_type<var>::type>();
}

TEST(MetaTraits, ReturnTypeFvarDouble) {
  test::expect_same_type<fvar<double>, return_type<fvar<double> >::type>();
}

TEST(MetaTraits, ReturnTypeFvarFvarDouble) {
  test::expect_same_type<fvar<fvar<double> >,
                         return_type<fvar<fvar<double> > >::type>();
}

TEST(MetaTraits, ReturnTypeFvarVar) {
  test::expect_same_type<fvar<var>, return_type<fvar<var> >::type>();
}

TEST(MetaTraits, ReturnTypeFvarFvarVar) {
  test::expect_same_type<fvar<fvar<var> >,
                         return_type<fvar<fvar<var> > >::type>();
}

TEST(MetaTraits, ReturnTypeScalarTenParams) {
  test::expect_same_type<double,
                         return_type<double, int, double, float, float, double,
                                     float, int, double, double>::type>();
}

TEST(MetaTraits, ReturnTypeVarTenParams) {
  test::expect_same_type<var,
                         return_type<double, var, double, int, double, float,
                                     float, float, var, int>::type>();
}

TEST(MetaTraits, ReturnTypeFvarDoubleTenParams) {
  test::expect_same_type<
      fvar<double>,
      return_type<double, fvar<double>, double, int, double, float, float,
                  float, fvar<double>, int>::type>();
}

TEST(MetaTraits, ReturnTypeFvarFvarDoubleTenParams) {
  test::expect_same_type<
      fvar<fvar<double> >,
      return_type<double, fvar<fvar<double> >, double, int, double, float,
                  float, float, fvar<fvar<double> >, int>::type>();
}

TEST(MetaTraits, ReturnTypeFvarVarTenParams) {
  test::expect_same_type<
      fvar<var>, return_type<double, fvar<var>, double, int, double, float,
                             float, float, fvar<var>, int>::type>();
}

TEST(MetaTraits, ReturnTypeFvarFvarVarTenParams) {
  test::expect_same_type<
      fvar<fvar<var> >,
      return_type<double, fvar<fvar<var> >, double, int, double, float, float,
                  float, fvar<fvar<var> >, int>::type>();
}
