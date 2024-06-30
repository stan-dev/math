#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, offset_multiplier) {
  EXPECT_FLOAT_EQ(2.0 - 5.0 * 1.0,
                  stan::math::offset_multiplier_constrain(-1.0, 2.0, 5.0));

  EXPECT_FLOAT_EQ(1.7, stan::math::offset_multiplier_constrain(1.7, 0, 1));
}

TEST(prob_transform, offset_multiplier_j) {
  double lp = -17.0;
  double L = 2.0;
  double U = 5.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(L + U * x,
                  stan::math::offset_multiplier_constrain(x, L, U, lp));
  EXPECT_FLOAT_EQ(-17.0 + log(U), lp);

  double lp1 = -12.9;
  EXPECT_FLOAT_EQ(1.7, stan::math::offset_multiplier_constrain(1.7, 0, 1, lp1));
  EXPECT_FLOAT_EQ(-12.9, lp1);
}

TEST(prob_transform, offset_multiplier_exception) {
  using stan::math::offset_multiplier_constrain;
  using stan::math::offset_multiplier_free;
  EXPECT_THROW(offset_multiplier_constrain(5.0, 1.0, 0.0), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(
                   5.0, std::numeric_limits<double>::infinity(), 1.0),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(5.0, NAN, 1.0), std::domain_error);
  EXPECT_NO_THROW(offset_multiplier_constrain(5.0, 1.0, 0.01));
  EXPECT_THROW(offset_multiplier_free(5.0, 1.0, 0.0), std::domain_error);
  EXPECT_THROW(
      offset_multiplier_free(5.0, std::numeric_limits<double>::infinity(), 1.0),
      std::domain_error);
  EXPECT_THROW(offset_multiplier_free(5.0, NAN, 1.0), std::domain_error);
  EXPECT_NO_THROW(offset_multiplier_free(5.0, 1.0, 0.01));
  double lp = 12;
  EXPECT_THROW(offset_multiplier_constrain(5.0, 1.0, 0.0, lp),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(
                   5.0, std::numeric_limits<double>::infinity(), 1.0, lp),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(5.0, NAN, 1.0, lp),
               std::domain_error);
  EXPECT_NO_THROW(offset_multiplier_constrain(5.0, 1.0, 0.01, lp));
}

TEST(prob_transform, offset_multiplier_f) {
  double L = -10.0;
  double U = 27.0;
  double y = 3.0;
  EXPECT_FLOAT_EQ(y, stan::math::offset_multiplier_constrain(
                         stan::math::offset_multiplier_free(y, L, U), L, U));
  EXPECT_FLOAT_EQ(y,
                  stan::math::offset_multiplier_free(
                      stan::math::offset_multiplier_constrain(y, L, U), L, U));
  L = 0.0;
  U = 1.0;
  y = 3.0;
  EXPECT_FLOAT_EQ(y, stan::math::offset_multiplier_constrain(
                         stan::math::offset_multiplier_free(y, L, U), L, U));
  EXPECT_FLOAT_EQ(y,
                  stan::math::offset_multiplier_free(
                      stan::math::offset_multiplier_constrain(y, L, U), L, U));
}

TEST(prob_transform, offset_multiplier_f_exception) {
  double L = -10.0;
  double U = -27.0;
  EXPECT_THROW(stan::math::offset_multiplier_free(L - 0.01, L, U),
               std::domain_error);
}

TEST(prob_transform, offset_multiplier_constrain_matrix) {
  Eigen::VectorXd x(4);
  x << -1.0, 1.1, 3.0, 5.0;
  Eigen::VectorXd sigma(4);
  sigma << 1.1, 0.3, 6.0, 7.0;
  Eigen::VectorXd offset(4);
  offset << -2.0, 0.0, 0.2, 2.0;

  double sigmad = 8.0;
  double offsetd = -2.0;

  Eigen::VectorXd sigma_bad(3);
  Eigen::VectorXd offset_bad(3);

  // matrix, real, real
  {
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offsetd, sigmad);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offsetd, sigmad));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offsetd, sigmad);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
  }
  // matrix, matrix, real
  {
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offset, sigmad);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offset(i), sigmad));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offset, sigmad);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offset_bad, sigmad),
                 std::invalid_argument);
  }
  // matrix, real, matrix
  {
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offsetd, sigma);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offsetd, sigma(i)));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offsetd, sigma);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offsetd, sigma_bad),
                 std::invalid_argument);
  }
  // matrix, matrix, matrix
  {
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offset, sigma);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offset(i), sigma(i)));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offset, sigma);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offset, sigma_bad),
                 std::invalid_argument);
  }

  // matrix, real, real, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offsetd, sigmad, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offsetd, sigmad, lp1));
    }
    EXPECT_EQ(lp0, lp1);
    auto x_free = stan::math::offset_multiplier_free(result, offsetd, sigmad);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
  }
  // matrix, matrix, real, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offset, sigmad, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offset(i), sigmad, lp1));
    }
    EXPECT_EQ(lp0, lp1);
    auto x_free = stan::math::offset_multiplier_free(result, offset, sigmad);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offset_bad, sigmad),
                 std::invalid_argument);
  }
  // matrix, real, matrix, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offsetd, sigma, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offsetd, sigma(i), lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    auto x_free = stan::math::offset_multiplier_free(result, offsetd, sigma);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offsetd, sigma_bad),
                 std::invalid_argument);
  }
  // matrix, matrix, matrix, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result
        = stan::math::offset_multiplier_constrain(x, offset, sigma, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::offset_multiplier_constrain(
                                     x(i), offset(i), sigma(i), lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    auto x_free = stan::math::offset_multiplier_free(result, offset, sigma);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offset, sigma_bad),
                 std::invalid_argument);
  }
}

TEST(prob_transform, offset_multiplier_constrain_stdvec) {
  std::vector<double> x{-1.0, 1.1, 3.0, 5.0};
  std::vector<double> sigma{1.1, 0.2, 6.0, 7.0};
  std::vector<double> offset{-2.0, 0.0, 0.2, 2.0};

  double sigmad = 8.0;
  double offsetd = -2.0;
  // These are for dimension checks
  std::vector<double> offset_bad{-2, -3};
  std::vector<double> sigma_bad{8, 9};

  // array[], real, real
  {
    auto result = stan::math::offset_multiplier_constrain(x, offsetd, sigmad);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::offset_multiplier_constrain(
                                     x[i], offsetd, sigmad));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offsetd, sigmad);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
  }
  // array[], array[], real
  {
    auto result = stan::math::offset_multiplier_constrain(x, offset, sigmad);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::offset_multiplier_constrain(
                                     x[i], offset[i], sigmad));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offset, sigmad);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offset_bad, sigmad),
                 std::invalid_argument);
  }
  // array[], real, array[]
  {
    auto result = stan::math::offset_multiplier_constrain(x, offsetd, sigma);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::offset_multiplier_constrain(
                                     x[i], offsetd, sigma[i]));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offsetd, sigma);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offsetd, sigma_bad),
                 std::invalid_argument);
  }
  // array[], array[], array[]
  {
    auto result = stan::math::offset_multiplier_constrain(x, offset, sigma);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::offset_multiplier_constrain(
                                     x[i], offset[i], sigma[i]));
    }
    auto x_free = stan::math::offset_multiplier_free(result, offset, sigma);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x, offset, sigma_bad),
                 std::invalid_argument);
  }
}

TEST(prob_transform, offset_multiplier_constrain_stdvec_matrix) {
  Eigen::VectorXd x(4);
  x << -1.0, 1.1, 3.0, 5.0;
  Eigen::VectorXd sigma(4);
  sigma << 1.1, 0.2, 6.0, 7.0;
  Eigen::VectorXd offset(4);
  offset << -2.0, 0.0, 0.2, 2.0;

  double sigmad = 8.0;
  double offsetd = -2.0;

  Eigen::VectorXd sigma_bad(3);
  Eigen::VectorXd offset_bad(3);

  std::vector<Eigen::VectorXd> x_vec{x, x};
  std::vector<Eigen::VectorXd> sigma_vec{sigma, sigma};
  std::vector<Eigen::VectorXd> offset_vec{offset, offset};
  std::vector<Eigen::VectorXd> sigma_bad_vec{sigma_bad, sigma_bad, sigma_bad};
  std::vector<Eigen::VectorXd> offset_bad_vec{offset_bad, offset_bad,
                                              sigma_bad};

  // array[] matrix, array[] matrix, array[] matrix
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offset_vec, sigma_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::offset_multiplier_constrain(
                            x_vec[i](j), offset_vec[i](j), sigma_vec[i](j)));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offset_vec, sigma_vec);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x_vec, offset_bad_vec,
                                                         sigma_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x_vec, offset_bad_vec,
                                                         sigma_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::offset_multiplier_constrain(x_vec, offset_vec,
                                                         sigma_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::offset_multiplier_free(x_vec, offset_bad_vec,
                                                    sigma_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offset_bad_vec, sigma_vec),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offset_vec, sigma_bad_vec),
        std::invalid_argument);
  }
  // array[] matrix, array[] matrix, real
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offset_vec, sigmad);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::offset_multiplier_constrain(
                            x_vec[i](j), offset_vec[i](j), sigmad));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offset_vec, sigmad);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offset_bad_vec, sigmad),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offset_bad_vec, sigmad),
        std::invalid_argument);
  }
  // array[] matrix, real, array[] matrix
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offsetd, sigma_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::offset_multiplier_constrain(
                            x_vec[i](j), offsetd, sigma_vec[i](j)));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offsetd, sigma_vec);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offsetd, sigma_bad_vec),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offsetd, sigma_bad_vec),
        std::invalid_argument);
  }
  // array[] matrix, array[] matrix, matrix
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offset_vec, sigma);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::offset_multiplier_constrain(
                            x_vec[i](j), offset_vec[i](j), sigma(j)));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offset_vec, sigma);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offset_vec, sigma_bad),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offset_vec, sigma_bad),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offset_bad_vec, sigma),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offset_bad_vec, sigma),
        std::invalid_argument);
  }
  // array[] matrix, matrix, array[] matrix
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offset, sigma_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::offset_multiplier_constrain(
                            x_vec[i](j), offset(j), sigma_vec[i](j)));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offset, sigma_vec);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offset_bad, sigma_vec),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offset_bad, sigma_vec),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offset, sigma_bad_vec),
        std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_free(x_vec, offset, sigma_bad_vec),
        std::invalid_argument);
  }
  // array[] matrix, matrix, matrix
  {
    auto result = stan::math::offset_multiplier_constrain(x_vec, offset, sigma);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::offset_multiplier_constrain(
                                          x_vec[i](j), offset(j), sigma(j)));
      }
    }
    auto x_vec_free = stan::math::offset_multiplier_free(result, offset, sigma);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offset_bad, sigma),
        std::invalid_argument);
    EXPECT_THROW(stan::math::offset_multiplier_free(x_vec, offset_bad, sigma),
                 std::invalid_argument);
    EXPECT_THROW(
        stan::math::offset_multiplier_constrain(x_vec, offset, sigma_bad),
        std::invalid_argument);
    EXPECT_THROW(stan::math::offset_multiplier_free(x_vec, offset, sigma_bad),
                 std::invalid_argument);
  }
  // The rest of the real versions of these don't need the throw test
  // array[] matrix, matrix, real
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offset, sigmad);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::offset_multiplier_constrain(
                                          x_vec[i](j), offset(j), sigmad));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offset, sigmad);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
  }
  // array[] matrix, real, matrix
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offsetd, sigma);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::offset_multiplier_constrain(
                                          x_vec[i](j), offsetd, sigma(j)));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offsetd, sigma);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
  }
  // array[] matrix, real, real
  {
    auto result
        = stan::math::offset_multiplier_constrain(x_vec, offsetd, sigmad);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::offset_multiplier_constrain(
                                          x_vec[i](j), offsetd, sigmad));
      }
    }
    auto x_vec_free
        = stan::math::offset_multiplier_free(result, offsetd, sigmad);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
  }
}

TEST(prob_transform, offset_multiplier_free_exception) {
  double xd = 1.0;
  Eigen::VectorXd x(2);
  x << 1.0, 2.0;

  double mu = 2.0;

  double mu_bad1 = stan::math::INFTY;
  Eigen::VectorXd mu_bad2(2);
  mu_bad2 << stan::math::INFTY, 2.0;
  Eigen::VectorXd mu_bad3(3);

  double sigma = 2.0;

  double sigma_bad1 = stan::math::INFTY;
  double sigma_bad2 = 0.0;
  double sigma_bad3 = -2.0;
  Eigen::VectorXd sigma_bad4(2);
  sigma_bad4 << stan::math::INFTY, 2.0;
  Eigen::VectorXd sigma_bad5(2);
  sigma_bad5 << -1.0, 2.0;
  Eigen::VectorXd sigma_bad6(3);

  EXPECT_THROW(stan::math::offset_multiplier_free(xd, mu_bad1, sigma),
               std::domain_error);

  EXPECT_THROW(stan::math::offset_multiplier_free(xd, mu, sigma_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(xd, mu, sigma_bad2),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(xd, mu, sigma_bad3),
               std::domain_error);

  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu, sigma_bad1),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu, sigma_bad2),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu, sigma_bad3),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu, sigma_bad4),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu, sigma_bad5),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu, sigma_bad6),
               std::invalid_argument);

  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu_bad1, sigma),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu_bad2, sigma),
               std::domain_error);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu_bad3, sigma),
               std::invalid_argument);
}

TEST(prob_transform, offset_multiplier_consistent_sizes) {
  double xd = 1.0;
  Eigen::VectorXd x(4);
  x << 1.0, 2.0, 3.0, 4.0;

  double mud = 2.0;
  Eigen::RowVectorXd mu(4);
  mu << 1.0, 2.0, 3.0, 4.0;

  double sigmad = 1.0;
  Eigen::MatrixXd sigma(2, 2);
  sigma << 1.0, 2.0, 3.0, 4.0;

  EXPECT_THROW(stan::math::offset_multiplier_constrain(xd, mu, sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::offset_multiplier_constrain(x, mud, sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::offset_multiplier_constrain(x, mu, sigmad),
               std::invalid_argument);

  EXPECT_THROW(stan::math::offset_multiplier_free(xd, mu, sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mud, sigma),
               std::invalid_argument);
  EXPECT_THROW(stan::math::offset_multiplier_free(x, mu, sigmad),
               std::invalid_argument);
}
