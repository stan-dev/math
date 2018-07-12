#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MetaTraits, ReturnTypeMix) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::return_type;

  test::expect_same_type<fvar<var>, return_type<fvar<var> >::type>();
  test::expect_same_type<fvar<var>, return_type<fvar<var>, fvar<var> >::type>();
  test::expect_same_type<fvar<var>,
                         return_type<fvar<var>, fvar<var>, fvar<var> >::type>();
  test::expect_same_type<fvar<var>,
                         return_type<fvar<var>, fvar<var>, var>::type>();

  test::expect_same_type<fvar<fvar<var> >,
                         return_type<fvar<fvar<var> > >::type>();
  test::expect_same_type<
      fvar<fvar<var> >,
      return_type<fvar<fvar<var> >, fvar<fvar<var> > >::type>();
  test::expect_same_type<fvar<fvar<var> >,
                         return_type<fvar<fvar<var> >, fvar<fvar<var> >,
                                     fvar<fvar<var> > >::type>();
  test::expect_same_type<
      fvar<fvar<var> >,
      return_type<fvar<fvar<var> >, fvar<fvar<var> >, var>::type>();
}
