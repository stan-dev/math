#include <stan/math/prim/scal.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::partials_type;

TEST(MetaTraits, PartialsTypeDouble) {
  test::expect_same_type<double, partials_type<double>::type>();
}

TEST(MetaTraits, PartialsTypeFloat) {
  test::expect_same_type<float, partials_type<float>::type>();
}

TEST(MetaTraits, PartialsTypeInt) {
  test::expect_same_type<int, partials_type<int>::type>();
}
