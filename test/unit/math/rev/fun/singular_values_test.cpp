#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/rev/fun/singular_values.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, singularvalues_gradient) {
  // logdet(A) can be calculated using singularvalues of matrix A
  // the derivative of logdet(A) should be inverse(A)

  Eigen::MatrixXd a(4, 4);
  // Random symmetric matrix
  a << 1.8904, 0.7204, -0.1599, 1.2028, 0.7204, 7.3394, 2.0895, -0.6151,
      -0.1599, 2.0895, 2.4601, -1.7219, 1.2028, -0.6151, -1.7219, 4.5260;
  stan::math::matrix_d a_inv = stan::math::inverse(a);

  stan::math::matrix_v a_v(a);
  auto w = stan::math::singular_values(a_v);
  auto logdet = stan::math::sum(stan::math::log(w));

  stan::math::set_zero_all_adjoints();
  logdet.grad();
  ASSERT_TRUE(a_inv.val().isApprox(a_v.adj()));
}
