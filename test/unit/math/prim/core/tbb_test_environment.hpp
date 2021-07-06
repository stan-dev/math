#ifdef TEST_UNIT_MATH_TBB_ENVIRONMENT
#define TEST_UNIT_MATH_TBB_ENVIRONMENT

#include <stan/math/prim/core/init_threadpool_tbb.hpp>

// Ensure that the Intel TBB is setup in alignment with STAN_THREADS.
class TBBEnvironment : public ::testing::Environment {
 public:
  // Initialise the Intel TBB threadpool
  virtual void SetUp() { stan::math::init_threadpool_tbb(); }
};

::testing::Environment* const tbb_env = ::testing::AddGlobalTestEnvironment(new TBBEnvironment);

#endif
