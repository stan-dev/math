#include <stan/math/mix/scal.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::var;
using stan::partials_return_type;

TEST(MetaTraits, PartialsReturnTypeFvarVar) {
  test::expect_same_type<var, partials_return_type<fvar<var> >::type>();
}

TEST(MetaTraits, PartialsReturnTypeFvarFvarVar) {
  test::expect_same_type<fvar<var>,
                         partials_return_type<fvar<fvar<var> > >::type>();
}

TEST(MetaTraits, PartialsReturnTypeFvarVarTenParams) {
  test::expect_same_type<
      var, partials_return_type<double, fvar<var>, double, int, double, float,
                                float, float, fvar<var>, int>::type>();
}

TEST(MetaTraits, PartialsReturnTypeFvarFvarVarTenParams) {
  test::expect_same_type<
      fvar<var>,
      partials_return_type<double, fvar<fvar<var> >, double, int, double, float,
                           float, float, fvar<fvar<var> >, int>::type>();
}
