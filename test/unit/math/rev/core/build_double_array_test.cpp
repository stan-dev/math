#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRev, vars_test) {
  std::vector<stan::math::var> vars = {1.0, 2.0, 3.0, -5.0, 10.0};

  stan::math::vari** varis = stan::math::build_vari_pointer_array_if_necessary(vars.data(), vars.size());
  
  const double* doubles
      = stan::math::build_double_array(varis, vars.size());

  EXPECT_TRUE(
      stan::math::ChainableStack::instance().memalloc_.in_stack(doubles));

  for (size_t i = 0; i < vars.size(); ++i)
    EXPECT_EQ(doubles[i], vars[i].val());
}

TEST(AgradRev, doubles_test) {
  std::vector<double> doubles_vector = {1.0, 2.0, 3.0, -5.0, 10.0};

  const double* doubles = stan::math::build_double_array(doubles_vector.data(), doubles_vector.size());

  EXPECT_NE(doubles, doubles_vector.data());

  for (size_t i = 0; i < doubles_vector.size(); ++i)
    EXPECT_EQ(doubles[i], doubles_vector[i]);
}
