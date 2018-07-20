#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRev, vars_test) {
  std::vector<stan::math::var> vars = {1.0, 2.0, 3.0, -5.0, 10.0};

  double* doubles
      = stan::math::build_double_array_if_necessary(vars.data(), vars.size());

  EXPECT_TRUE(
      stan::math::ChainableStack::instance().memalloc_.in_stack(doubles));

  for (int i = 0; i < vars.size(); ++i)
    EXPECT_EQ(doubles[i], vars[i].val());
}

TEST(AgradRev, doubles_test) {
  std::vector<double> doubles_vector = {1.0, 2.0, 3.0, -5.0, 10.0};

  double* doubles = stan::math::build_double_array_if_necessary(
      doubles_vector.data(), doubles_vector.size());

  EXPECT_EQ(doubles, doubles_vector.data());

  for (int i = 0; i < doubles_vector.size(); ++i)
    EXPECT_EQ(doubles[i], doubles_vector[i]);
}
