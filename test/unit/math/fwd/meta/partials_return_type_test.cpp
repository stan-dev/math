#include <stan/math/fwd/arr.hpp>
#include <stan/math/fwd/scal.hpp>
#include <stan/math/fwd/meta.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

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

TEST(MetaTraitsFwdArr, partials_return_type) {
  partials_return_type<double, fvar<double>, std::vector<fvar<double> > >::type
      a(5.0);
  EXPECT_EQ(5.0, a);
}
