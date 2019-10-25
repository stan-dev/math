#include <stan/math/fwd/scal.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::partials_return_type;

TEST(MetaTraitsFwdScal, PartialsReturnTypeFvarDouble) {
  test::expect_same_type<double, partials_return_type<fvar<double> >::type>();
}

TEST(MetaTraitsFwdScal, PartialsReturnTypeFvarFvarDouble) {
  test::expect_same_type<fvar<double>,
                         partials_return_type<fvar<fvar<double> > >::type>();
}

TEST(MetaTraitsFwdScal, PartialsReturnTypeFvarDoubleTenParams) {
  test::expect_same_type<
      double,
      partials_return_type<double, fvar<double>, double, int, double, float,
                           float, float, fvar<double>, int>::type>();
}

TEST(MetaTraitsFwdScal, PartialsReturnTypeFvarFvarDoubleTenParams) {
  test::expect_same_type<
      fvar<double>, partials_return_type<double, fvar<fvar<double> >, double,
                                         int, double, float, float, float,
                                         fvar<fvar<double> >, int>::type>();
}
