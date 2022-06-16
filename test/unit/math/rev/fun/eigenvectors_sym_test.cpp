#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/rev/fun/eigenvectors_sym.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, eigenvectorsSym) {
  Eigen::MatrixXd a(4, 4);
  // Random symmetric matrix
  a << 1.8904, 0.7204, -0.1599, 1.2028, 0.7204, 7.3394, 2.0895, -0.6151,
      -0.1599, 2.0895, 2.4601, -1.7219, 1.2028, -0.6151, -1.7219, 4.5260;

  stan::math::matrix_v a_w(a);
  auto w = eigenvectors_sym(a_w);

  // TODO(Simon): how to test this?
}
