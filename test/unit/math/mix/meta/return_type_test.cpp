#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::var;
using stan::return_type;

TEST(MetaTraitsMixScal, ReturnTypeFvarVar) {
  test::expect_same_type<fvar<var>, return_type<fvar<var> >::type>();
}

TEST(MetaTraitsMixScal, ReturnTypeFvarFvarVar) {
  test::expect_same_type<fvar<fvar<var> >,
                         return_type<fvar<fvar<var> > >::type>();
}

TEST(MetaTraitsMixScal, ReturnTypeFvarVarTenParams) {
  test::expect_same_type<
      fvar<var>, return_type<double, fvar<var>, double, int, double, float,
                             float, float, fvar<var>, int>::type>();
}

TEST(MetaTraitsMixScal, ReturnTypeFvarFvarVarTenParams) {
  test::expect_same_type<
      fvar<fvar<var> >,
      return_type<double, fvar<fvar<var> >, double, int, double, float, float,
                  float, fvar<fvar<var> >, int>::type>();
}
