#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, lb_constrain) {
  double x = -1.0;
  double lb = 2.0;
  double lbi = stan::math::NEGATIVE_INFTY;

  {
    double result = stan::math::lb_constrain(x, lb);
    double resulti = stan::math::lb_constrain(x, lbi);
    EXPECT_FLOAT_EQ(exp(x) + lb, result);
    EXPECT_FLOAT_EQ(x, resulti);
    EXPECT_FLOAT_EQ(x, stan::math::lb_free(result, lb));
    EXPECT_FLOAT_EQ(x, stan::math::lb_free(resulti, lbi));
  }

  {
    double lp = 1.0;
    double lpi = 1.0;
    double result = stan::math::lb_constrain(x, lb, lp);
    double resulti = stan::math::lb_constrain(x, lbi, lpi);
    EXPECT_FLOAT_EQ(exp(x) + lb, result);
    EXPECT_FLOAT_EQ(x, resulti);
    EXPECT_FLOAT_EQ(0.0, lp);
    EXPECT_FLOAT_EQ(1.0, lpi);
  }
}

TEST(prob_transform, lb_constrain_matrix) {
  Eigen::VectorXd x(2);
  x << -1.0, 1.1;
  Eigen::VectorXd lb(2);
  lb << 2.0, stan::math::NEGATIVE_INFTY;
  double lbd = 2.0;

  Eigen::VectorXd lb_bad(3);

  // matrix, real
  {
    Eigen::VectorXd result = stan::math::lb_constrain(x, lbd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lb_constrain(x(i), lbd));
    }
    auto x_free = stan::math::lb_free(result, lbd);
    // EXPECT_MATRIX_EQ(x, x_free);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
  }
  // matrix, matrix
  {
    Eigen::VectorXd result = stan::math::lb_constrain(x, lb);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lb_constrain(x(i), lb(i)));
    }
    EXPECT_MATRIX_EQ(x, stan::math::lb_free(result, lb));
    EXPECT_THROW(stan::math::lb_constrain(x, lb_bad), std::invalid_argument);
  }
  // matrix, real
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::lb_constrain(x, lbd, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lb_constrain(x(i), lbd, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
  }
  // matrix, matrix
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::lb_constrain(x, lb, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lb_constrain(x(i), lb(i), lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    EXPECT_THROW(stan::math::lb_constrain(x, lb_bad, lp1),
                 std::invalid_argument);
  }
}

TEST(prob_transform, lb_constrain_std_vector) {
  Eigen::VectorXd x(2);
  x << -1.0, 1.1;
  Eigen::VectorXd lb(2);
  lb << 2.0, stan::math::NEGATIVE_INFTY;
  Eigen::VectorXd lb2(2);
  lb << stan::math::NEGATIVE_INFTY, 1.0;
  double lbd = 2.0;

  Eigen::VectorXd lb_bad(3);

  std::vector<Eigen::VectorXd> x_vec = {x, x};
  std::vector<Eigen::VectorXd> lb_vec = {lb, lb2};

  std::vector<Eigen::VectorXd> lb_vec_bad1 = {lb, lb2, lb};
  std::vector<Eigen::VectorXd> lb_vec_bad2 = {lb, lb_bad};
  // matrix[], real
  {
    std::vector<Eigen::VectorXd> result = stan::math::lb_constrain(x_vec, lbd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::lb_constrain(x_vec[i], lbd));
    }
    auto free_x = stan::math::lb_free(result, lbd);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        EXPECT_FLOAT_EQ(x.coeff(j), free_x[i].coeff(j));
      }
    }
  }

  // matrix[], matrix
  {
    std::vector<Eigen::VectorXd> result = stan::math::lb_constrain(x_vec, lb);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::lb_constrain(x_vec[i], lb));
    }
    auto free_x = stan::math::lb_free(result, lb);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        EXPECT_FLOAT_EQ(x.coeff(j), free_x[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lb_constrain(x_vec, lb_bad),
                 std::invalid_argument);
  }

  // matrix[], matrix[]
  {
    std::vector<Eigen::VectorXd> result
        = stan::math::lb_constrain(x_vec, lb_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i],
                       stan::math::lb_constrain(x_vec[i], lb_vec[i]));
    }
    auto free_x = stan::math::lb_free(result, lb_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < x.size(); ++j) {
        EXPECT_FLOAT_EQ(x.coeff(j), free_x[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lb_constrain(x_vec, lb_vec_bad1),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lb_constrain(x_vec, lb_vec_bad2),
                 std::invalid_argument);
  }

  // matrix[], real
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    std::vector<Eigen::VectorXd> result
        = stan::math::lb_constrain(x_vec, lbd, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::lb_constrain(x_vec[i], lbd, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
  }

  // matrix[], matrix
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    std::vector<Eigen::VectorXd> result
        = stan::math::lb_constrain(x_vec, lb, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i], stan::math::lb_constrain(x_vec[i], lb, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    EXPECT_THROW(stan::math::lb_constrain(x_vec, lb_bad, lp1),
                 std::invalid_argument);
  }

  // matrix[], matrix[]
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    std::vector<Eigen::VectorXd> result
        = stan::math::lb_constrain(x_vec, lb_vec, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_MATRIX_EQ(result[i],
                       stan::math::lb_constrain(x_vec[i], lb_vec[i], lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    EXPECT_THROW(stan::math::lb_constrain(x_vec, lb_vec_bad1, lp1),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lb_constrain(x_vec, lb_vec_bad2, lp1),
                 std::invalid_argument);
  }
}

TEST(prob_transform, lb_free_exception) {
  double lbd = 2.0;
  Eigen::VectorXd lb(2);
  lb << 1.0, 2.0;
  Eigen::VectorXd lb2(2);
  lb2 << -1.0, 5.0;
  std::vector<Eigen::VectorXd> lb_vec = {lb, lb2};

  double xd_bad = -5.0;
  Eigen::VectorXd x_bad(2);
  x_bad << -10.0, -5.0;
  Eigen::VectorXd x_bad2(2);
  x_bad2 << -10.0, -5.0;
  std::vector<Eigen::VectorXd> x_bad_vec = {x_bad, x_bad2};

  EXPECT_THROW(stan::math::lb_free(xd_bad, lbd), std::domain_error);

  EXPECT_THROW(stan::math::lb_free(x_bad, lbd), std::domain_error);
  EXPECT_THROW(stan::math::lb_free(x_bad, lb), std::domain_error);

  EXPECT_THROW(stan::math::lb_free(x_bad_vec, lbd), std::domain_error);
  EXPECT_THROW(stan::math::lb_free(x_bad_vec, lb), std::domain_error);
  EXPECT_THROW(stan::math::lb_free(x_bad_vec, lb_vec), std::domain_error);
}
