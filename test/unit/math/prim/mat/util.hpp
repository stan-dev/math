#ifndef TEST_MATH_UNIT_PRIM_MAT_UTIL_HPP
#define TEST_MATH_UNIT_PRIM_MAT_UTIL_HPP

#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

namespace stan {
  namespace test {
    namespace unit {

      void expect_symmetric(const Eigen::MatrixXd& a) {
        for (int j = 1; j < a.cols(); ++j)
          for (int i = 0; i < j; ++i)
            EXPECT_EQ(a(i, j), a(j, i))
              << "failed symmetry at " << i << ", " << j;
      }

    }
  }
}
#endif
