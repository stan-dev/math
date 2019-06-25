
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>






using stan::partials_type;

TEST(MetaTraits_scal, PartialsTypeDouble) {
  test::expect_same_type<double, partials_type<double>::type>();
}

TEST(MetaTraits_scal, PartialsTypeFloat) {
  test::expect_same_type<float, partials_type<float>::type>();
}

TEST(MetaTraits_scal, PartialsTypeInt) {
  test::expect_same_type<int, partials_type<int>::type>();
}
