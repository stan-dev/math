#include <stan/math/prim/arr.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, ScalarTypeArray) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double, scalar_type<vector<double>>::type>();
  test::expect_same_type<int, scalar_type<vector<int>>::type>();
  test::expect_same_type<double, scalar_type<vector<vector<double>>>::type>();
}

TEST(MetaTraits, ScalarTypeArrayConst) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double, scalar_type<const vector<double>>::type>();
  test::expect_same_type<int, scalar_type<const vector<int>>::type>();
  test::expect_same_type<double,
                         scalar_type<const vector<vector<double>>>::type>();
}

TEST(MetaTraits, ScalarTypeArrayConstConst) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double const*,
                         scalar_type<const vector<double const*>>::type>();
  test::expect_same_type<int const*,
                         scalar_type<const vector<int const*>>::type>();
  test::expect_same_type<
      double const*, scalar_type<const vector<vector<double const*>>>::type>();
}
