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
                   std::domain_error,
                   "An invalid Butcher or MRI table was encountered");
}

TEST(ErrorHandling, cvodes_check) {
  EXPECT_NO_THROW(stan::math::cvodes_check(CV_SUCCESS, "CVode"));
  EXPECT_THROW_MSG(stan::math::cvodes_check(CV_TOO_MUCH_WORK, "CVode"),
                   std::domain_error, "CVode failed with error flag -1");
  EXPECT_THROW_MSG(stan::math::cvodes_check(CV_BAD_IS, "f"), std::domain_error,
                   "The sensitivity index is larger than the number of "
                   "sensitivities computed");
  EXPECT_EQ(stan::math::cvodes_flag_msg(CV_NO_SENS)[0], "CV_NO_SENS");
  EXPECT_EQ(stan::math::cvodes_flag_msg(-1000)[0], "CV_UNKNOWN_ERROR");
}

TEST(ErrorHandling, idas_check) {
  EXPECT_NO_THROW(stan::math::idas_check(IDA_SUCCESS, "IDASolve"));
  EXPECT_THROW_MSG(stan::math::idas_check(IDA_TOO_MUCH_WORK, "IDASolve"),
                   std::domain_error, "IDASolve failed with error flag -1");
  EXPECT_THROW_MSG(stan::math::idas_check(IDA_BAD_EWT, "f"), std::domain_error,
                   "Zero value of some error weight component.");
  EXPECT_EQ(stan::math::idas_flag_msg(IDA_REP_RES_ERR)[0], "IDA_REP_RES_ERR");
  EXPECT_EQ(stan::math::idas_flag_msg(-1000)[0], "IDA_UNKNOWN_ERROR");
}
