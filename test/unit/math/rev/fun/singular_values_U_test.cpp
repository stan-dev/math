#include <stan/math/prim.hpp>
#include <stan/math/rev/fun/singular_values_U.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <typeinfo>
TEST(AgradRev, svd_U_gradient) {
    using namespace stan::math;
    Eigen::MatrixXd A(4, 4);
    // Random symmetric matrix
    A << 1.8904, 0.7204, -0.1599, 1.2028, 0.7204, 7.3394, 2.0895, -0.6151,
    -0.1599, 2.0895, 2.4601, -1.7219, 1.2028, -0.6151, -1.7219, 4.5260;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    auto U = svd.matrixU();
    Eigen::MatrixXd D = svd.singularValues().asDiagonal();
    auto V = svd.matrixV();

    Eigen:: MatrixXd Fp(4, 4);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i == j){
                Fp(i, j) = 0.0;
            }
            else{
                Fp(i, j) = 1.0 / (D(j, j) - D(i, i)) + 1.0 / (D(i, i) + D(j, j));
            }
        }
    }

    Eigen::MatrixXd adj_A_fd = Eigen::MatrixXd::Zero(4, 4);

    Eigen::MatrixXd adj_U_fd = Eigen::MatrixXd::Ones(4, 4);

    adj_A_fd += 0.5 * U * 
        (Fp.array() * (U.transpose() * adj_A_fd - adj_A_fd.transpose() * U).array()).matrix()
        + (Eigen::MatrixXd::Identity(4, 4) - U * U.transpose()) * adj_U_fd * D.inverse() * V.transpose();


    stan::math::set_zero_all_adjoints();
    matrix_v arena_a(A);
    auto stan_u = svd_U(arena_a);
    EXPECT_MATRIX_FLOAT_EQ(U, stan_u.val());
    EXPECT_MATRIX_FLOAT_EQ(adj_A_fd, arena_a.adj_op());
}