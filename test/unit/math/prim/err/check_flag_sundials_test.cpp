#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandling, arkode_check) {
  EXPECT_NO_THROW(stan::math::arkode_check(ARK_SUCCESS, "ERKStepEvolve"));
  EXPECT_NO_THROW(stan::math::arkode_check(ARK_TSTOP_RETURN, "ERKStepEvolve"));
  EXPECT_THROW_MSG(stan::math::arkode_check(ARK_TOO_MUCH_ACC, "ERKStepEvolve"),
                   std::domain_error,
                   "ERKStepEvolve failed with error flag -2");
  EXPECT_THROW_MSG(stan::math::arkode_check(ARK_INVALID_TABLE, "f"),
                   std::domain_error, "The Butcher table is invalid");
}
