#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <vector>

TEST(MathPrimMat, vec_double_dot_product0) {
    std::vector<double> x0(3);
    x0[0] = 13;
    x0[1] = -17;
    x0[2] = 19;

    std::vector<double> x1(3);
    x1[0] = -3;
    x1[1] = 1;
    x1[2] = 9;

    double dot_product;
    dot_product = stan::math::dot_product(x0, x1, 3);

    double answer = 0;
    for (int i = 0; i < 3; ++i) {
        answer += x0[i] * x1[i];
    }
    EXPECT_FLOAT_EQ(answer, dot_product);
}

TEST(MathPrimMat, vec_double_float_dot_product0) {
    std::vector<double> x0(3);
    x0[0] = 13.0;
    x0[1] = -17.0;
    x0[2] = 19.0;

    std::vector<float> x1(3);
    x1[0] = -3.0;
    x1[1] = 1.0;
    x1[2] = 9.0;

    double dot_product;
    dot_product = stan::math::dot_product(x0, x1, 3);

    double answer = 0;
    for (int i = 0; i < 3; ++i) {
        answer += x0[i] * x1[i];
    }
    EXPECT_FLOAT_EQ(answer, dot_product);
}

TEST(MathPrimMat, vec_int_float_dot_product0) {
    std::vector<int> x0(3);
    x0[0] = 13;
    x0[1] = -17;
    x0[2] = 19;

    std::vector<float> x1(3);
    x1[0] = -3.0;
    x1[1] = 1.0;
    x1[2] = 9.0;

    double dot_product;
    dot_product = stan::math::dot_product(x0, x1, 3);

    double answer = 0;
    for (int i = 0; i < 3; ++i) {
        answer += x0[i] * x1[i];
    }
    EXPECT_FLOAT_EQ(answer, dot_product);
}

TEST(MathPrimMat, matrix_double_dot_product1) {
    std::vector<Eigen::Matrix<double, 1, -1>> x1(3);
    std::vector<Eigen::Matrix<double, 1, -1>> x2(3);
    double answer0;
    double answer1;
    double answer2;

    for (size_t i = 0; i < x1.size(); ++i) {
        x1[i].resize(1, 3);
        x1[i] << 1 * i, 2 * i, 3 * i;
    }

    for (size_t i = 0; i < x2.size(); ++i) {
        x2[i].resize(1, 3);
        x2[i] << 2 * i, 3 * i, 4 * i;
    }

    EXPECT_NO_THROW(answer0 = stan::math::dot_product(x1[0], x2[0]));
    EXPECT_NO_THROW(answer1 = stan::math::dot_product(x1[1], x2[1]));
    EXPECT_NO_THROW(answer2 = stan::math::dot_product(x1[2], x2[2]));

    ASSERT_FLOAT_EQ(0, answer0);
    ASSERT_FLOAT_EQ(20, answer1);
    ASSERT_FLOAT_EQ(80, answer2);
}

TEST(MathPrimMat, dot_prod_same_function_squared_distance) {
    std::vector<double> x(3);
    x[0] = -2;
    x[1] = -1;
    x[2] = -0.5;

    EXPECT_NO_THROW(stan::math::squared_distance(x[0], x[1]));
    EXPECT_NO_THROW(stan::math::dot_product(x[0], x[1]));
    EXPECT_NO_THROW(stan::math::dot_product(x, x));

    std::vector<Eigen::Matrix<double, -1, 1>> x1(3);
    for (size_t i = 0; i < x1.size(); ++i) {
        x1[i].resize(3, 1);
        x1[i] << 1 * i, 2 * i, 3 * i;
    }
    std::cout << x1.size() << "\n";
    std::cout << x[0].size() << "\n";

    EXPECT_NO_THROW(stan::math::squared_distance(x1[0], x1[1]));
    EXPECT_NO_THROW(stan::math::dot_product(x1[0], x1[1]));
    EXPECT_NO_THROW(stan::math::cov_exp_quad(x, 1, 1));
    EXPECT_NO_THROW(stan::math::squared_distance(1, 2));
}
