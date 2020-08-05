#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, AD_stack_matrix_matrix_test) {
    using stan::math::AD_stack_matrix;
    using Eigen::MatrixXd;

    //construction
    AD_stack_matrix<MatrixXd> a;
    AD_stack_matrix<MatrixXd> a2;
    AD_stack_matrix<MatrixXd> b(3,2);
    AD_stack_matrix<MatrixXd> b2(3,2);
    AD_stack_matrix<MatrixXd> c(MatrixXd::Ones(3,2));
    AD_stack_matrix<MatrixXd> d(c);
    AD_stack_matrix<MatrixXd> e(2*d);

    //assignment
    a = c;
    a2 = std::move(d);
    b = d;
    b2 = std::move(c);
    e = e + a;
    a = MatrixXd::Ones(3,2);

    EXPECT_MATRIX_EQ(a + a2 + b + b2 + e, MatrixXd::Ones(3,2) * 7);
}

TEST(AgradRev, AD_stack_matrix_vector_test) {
    using stan::math::AD_stack_matrix;
    using Eigen::VectorXd;

    //construction
    AD_stack_matrix<VectorXd> a;
    AD_stack_matrix<VectorXd> a2;
    AD_stack_matrix<VectorXd> b(3);
    AD_stack_matrix<VectorXd> b2(3);
    AD_stack_matrix<VectorXd> c(VectorXd::Ones(3));
    AD_stack_matrix<VectorXd> d(c);
    AD_stack_matrix<VectorXd> e(2*d);

    //assignment
    a = c;
    a2 = std::move(d);
    b = d;
    b2 = std::move(c);
    e = e + a;
    a = VectorXd::Ones(3);

    EXPECT_MATRIX_EQ(a + a2 + b + b2 + e, VectorXd::Ones(3) * 7);
}

TEST(AgradRev, AD_stack_matrix_row_vector_test) {
    using stan::math::AD_stack_matrix;
    using Eigen::RowVectorXd;

    //construction
    AD_stack_matrix<RowVectorXd> a;
    AD_stack_matrix<RowVectorXd> a2;
    AD_stack_matrix<RowVectorXd> b(3);
    AD_stack_matrix<RowVectorXd> b2(3);
    AD_stack_matrix<RowVectorXd> c(RowVectorXd::Ones(3));
    AD_stack_matrix<RowVectorXd> d(c);
    AD_stack_matrix<RowVectorXd> e(2*d);

    //assignment
    a = c;
    a2 = std::move(d);
    b = d;
    b2 = std::move(c);
    e = e + a;
    a = RowVectorXd::Ones(3);

    EXPECT_MATRIX_EQ(a + a2 + b + b2 + e, RowVectorXd::Ones(3) * 7);
}
