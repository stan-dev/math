#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, ub_constrain) {
  double x = -1.0;
  double ub = 2.0;
  double ubi = stan::math::INFTY;

  {
    double result = stan::math::ub_constrain(x, ub);
    double resulti = stan::math::ub_constrain(x, ubi);
    EXPECT_FLOAT_EQ(ub - exp(x), result);
    EXPECT_FLOAT_EQ(x, resulti);
    EXPECT_FLOAT_EQ(x, stan::math::ub_free(result, ub));
    EXPECT_FLOAT_EQ(x, stan::math::ub_free(resulti, ubi));
  }

  {
    double lp = 1.0;
    double lpi = 1.0;
    double result = stan::math::ub_constrain(x, ub, lp);
    double resulti = stan::math::ub_constrain(x, ubi, lpi);
    EXPECT_FLOAT_EQ(ub - exp(x), result);
    EXPECT_FLOAT_EQ(x, resulti);
    EXPECT_FLOAT_EQ(0.0, lp);
    EXPECT_FLOAT_EQ(1.0, lpi);
  }
}

TEST(prob_transform, ub_constrain_matrix) {
  Eigen::VectorXd x(2);
  x << -1.0, 1.1;
  Eigen::VectorXd ub(2);
  ub << 2.0, stan::math::INFTY;
  double ubd = 2.0;

  Eigen::VectorXd ub_bad(3);

  // matrix, real
  {
    Eigen::VectorXd result = stan::math::ub_constrain(x, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::ub_constrain(x(i), ubd));
    }
    auto x_free = stan::math::ub_free(result, ubd);
    // EXPECT_MATRIX_EQ(x, x_free);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
  }
  // matrix, matrix
  {
    Eigen::VectorXd result = stan::math::ub_constrain(x, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::ub_constrain(x(i), ub(i)));
    }
    EXPECT_MATRIX_EQ(x, stan::math::ub_free(result, ub));
    EXPECT_THROW(stan::math::ub_constrain(x, ub_bad), std::invalid_argument);
  }
  // matrix, real
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::ub_constrain(x, ubd, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::ub_constrain(x(i), ubd, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
  }

  // matrix, matrix
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::ub_constrain(x, ub, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::ub_constrain(x(i), ub(i), lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    EXPECT_THROW(stan::math::ub_constrain(x, ub_bad, lp1),
                 std::invalid_argument);
  }
}

TEST(prob_transform, ub_constrain_std_vector) {
  Eigen::VectorXd x(2);
  x << -1.0, 1.1;
  Eigen::VectorXd ub(2);
  ub << 2.0, stan::math::INFTY;
  Eigen::VectorXd ub2(2);
  ub2 << stan::math::INFTY, 1.0;
  double ubd = 2.0;

  Eigen::VectorXd ub_bad(3);

  std::vector<Eigen::VectorXd> x_vec = {x, 2 * x};
  std::vector<Eigen::VectorXd> ub_vec = {ub, ub2};

  std::vector<Eigen::VectorXd> ub_vec_bad1 = {ub, ub2, ub};
  std::vector<Eigen::VectorXd> ub_vec_bad2 = {ub, ub_bad};

  // matrix[], real
  {
    std::vector<Eigen::VectorXd> result = stan::math::ub_constrain(x_vec, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::ub_constrain(x_vec[i], ubd));
    }
    auto free_x = stan::math::ub_free(result, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), free_x[i].coeff(j));
      }
    }
  }

  // matrix[], matrix
  {
    std::vector<Eigen::VectorXd> result = stan::math::ub_constrain(x_vec, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::ub_constrain(x_vec[i], ub));
    }
    auto free_x = stan::math::ub_free(result, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), free_x[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::ub_constrain(x_vec, ub_bad),
                 std::invalid_argument);
  }

  // matrix[], matrix[]
  {
    std::vector<Eigen::VectorXd> result
        = stan::math::ub_constrain(x_vec, ub_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i],
                       stan::math::ub_constrain(x_vec[i], ub_vec[i]));
    }
    auto free_x = stan::math::ub_free(result, ub_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), free_x[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::ub_constrain(x_vec, ub_vec_bad1),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::ub_constrain(x_vec, ub_vec_bad2),
                 std::invalid_argument);
  }

  // matrix[], real
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    std::vector<Eigen::VectorXd> result
        = stan::math::ub_constrain(x_vec, ubd, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::ub_constrain(x_vec[i], ubd, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
  }

  // matrix[], matrix
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    std::vector<Eigen::VectorXd> result
        = stan::math::ub_constrain(x_vec, ub, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::ub_constrain(x_vec[i], ub, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    EXPECT_THROW(stan::math::ub_constrain(x_vec, ub_bad, lp1),
                 std::invalid_argument);
  }

  // matrix[], matrix[]
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    std::vector<Eigen::VectorXd> result
        = stan::math::ub_constrain(x_vec, ub_vec, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i],
                       stan::math::ub_constrain(x_vec[i], ub_vec[i], lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    EXPECT_THROW(stan::math::ub_constrain(x_vec, ub_vec_bad1, lp1),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::ub_constrain(x_vec, ub_vec_bad2, lp1),
                 std::invalid_argument);
  }
}

TEST(prob_transform, ub_free_exception) {
  double ubd = 2.0;
  Eigen::VectorXd ub(2);
  ub << 1.0, 2.0;
  Eigen::VectorXd ub2(2);
  ub2 << -1.0, 5.0;
  std::vector<Eigen::VectorXd> ub_vec = {ub, ub2};

  double xd_bad = 5.0;
  Eigen::VectorXd x_bad(2);
  x_bad << 10.0, 5.0;
  Eigen::VectorXd x_bad2(2);
  x_bad2 << 10.0, 5.0;
  std::vector<Eigen::VectorXd> x_bad_vec = {x_bad, x_bad2};

  EXPECT_THROW(stan::math::ub_free(xd_bad, ubd), std::domain_error);

  EXPECT_THROW(stan::math::ub_free(x_bad, ubd), std::domain_error);
  EXPECT_THROW(stan::math::ub_free(x_bad, ub), std::domain_error);

  EXPECT_THROW(stan::math::ub_free(x_bad_vec, ubd), std::domain_error);
  EXPECT_THROW(stan::math::ub_free(x_bad_vec, ub), std::domain_error);
  EXPECT_THROW(stan::math::ub_free(x_bad_vec, ub_vec), std::domain_error);
}
