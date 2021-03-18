#ifndef STAN_TEST_UNIT_MATH_PRIM_ERR_UTIL_HPP
#define STAN_TEST_UNIT_MATH_PRIM_ERR_UTIL_HPP

#ifdef STAN_NO_RANGE_CHECKS
/**
 * When STAN_NO_RANGE_CHECKS is defined this causes exceptions to turn off
 */
#define STAN_EXPECT_THROW(statement, err_type) EXPECT_NO_THROW(statement)
#define STAN_EXPECT_THROW_MSG(statement, err_type, message) \
  EXPECT_NO_THROW(statement)
#else
#define STAN_EXPECT_THROW(statement, err_type) EXPECT_THROW(statement, err_type)
#define STAN_EXPECT_THROW_MSG(statement, err_type, message) \
  EXPECT_THROW_MSG(statement, err_type, message)
#endif

#endif
