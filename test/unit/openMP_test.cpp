#include <gtest/gtest.h>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>

TEST(openMP, openMP_loop) {
  omp_set_num_threads(2);  // ordinarily picked up from environmental variable
  EXPECT_EQ(1, omp_get_num_threads());  // now outside parallel zone
  bool runtime = true;  // ordinarily a genuine runtime condition like N > 100
// We put pragmas like this before loops only if autodiff is not used
#pragma omp parallel for default(none) shared(std::cout) if (runtime)
  for (int i = 0; i < 2; i++) {
    EXPECT_EQ(2, omp_get_num_threads());  // now inside parallel zone
    std::cout << "Hello from thread " << omp_get_thread_num() << "\n"
              << std::endl;
  }
}
#endif
