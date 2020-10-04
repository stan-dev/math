#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>

TEST(AgradRevMatrix, multiply_lower_tri_self_transpose_test_large) {
  for (auto N : {1, 2, 3, 5, 10}) {
    Eigen::MatrixXd A
        = Eigen::MatrixXd::Random(N, N).template triangularView<Eigen::Lower>();
    stan::math::promote_scalar_t<stan::math::var, Eigen::MatrixXd> Av = A;

    auto res = multiply_lower_tri_self_transpose(Av);

    stan::math::sum(res).grad();

    Eigen::MatrixXd adjoint = (2 * Eigen::MatrixXd::Ones(N, N) * A)
                                  .template triangularView<Eigen::Lower>();
    EXPECT_MATRIX_NEAR(res.val(), A * A.transpose(), 1e-8);
    EXPECT_MATRIX_NEAR(Av.adj(), adjoint, 1e-8);

    stan::math::recover_memory();
  }
}
