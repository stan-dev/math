#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, lub) {
  EXPECT_FLOAT_EQ(2.0 + (5.0 - 2.0) * stan::math::inv_logit(-1.0),
                  stan::math::lub_constrain(-1.0, 2.0, 5.0));
}

TEST(prob_transform, lub_vec) {
  Eigen::VectorXd input(2);
  input << -1.0, 1.1;
  Eigen::VectorXd lbv(2);
  lbv << -1.0, 0.5;
  Eigen::VectorXd ubv(2);
  ubv << 2.0, 3.0;
  double lb = 1.0;
  double ub = 2.0;

  Eigen::VectorXd resvv(2);
  resvv << (2.0 + 1.0) * stan::math::inv_logit(-1.0) - 1.0,
      (3.0 - 0.5) * stan::math::inv_logit(1.1) + 0.5;
  Eigen::VectorXd ressv(2);
  ressv << (2.0 - 1.0) * stan::math::inv_logit(-1.0) + 1.0,
      (3.0 - 1.0) * stan::math::inv_logit(1.1) + 1.0;
  Eigen::VectorXd resvs(2);
  resvs << (2.0 + 1.0) * stan::math::inv_logit(-1.0) - 1.0,
      (2.0 - 0.5) * stan::math::inv_logit(1.1) + 0.5;
  Eigen::VectorXd res(2);
  res << (2.0 - 1.0) * stan::math::inv_logit(-1.0) + 1.0,
      (2.0 - 1.0) * stan::math::inv_logit(1.1) + 1.0;

  EXPECT_MATRIX_EQ(resvv, stan::math::lub_constrain(input, lbv, ubv));
  EXPECT_MATRIX_EQ(ressv, stan::math::lub_constrain(input, lb, ubv));
  EXPECT_MATRIX_EQ(resvs, stan::math::lub_constrain(input, lbv, ub));
  EXPECT_MATRIX_EQ(res, stan::math::lub_constrain(input, lb, ub));

  double lp = 0.0;
  EXPECT_MATRIX_EQ(resvv, stan::math::lub_constrain(input, lbv, ubv, lp));
  EXPECT_FLOAT_EQ(((ubv - lbv).array().log() + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
  lp = 0.0;
  EXPECT_MATRIX_EQ(ressv, stan::math::lub_constrain(input, lb, ubv, lp));
  EXPECT_FLOAT_EQ(((ubv.array() - lb).log() + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
  lp = 0.0;
  EXPECT_MATRIX_EQ(resvs, stan::math::lub_constrain(input, lbv, ub, lp));
  EXPECT_FLOAT_EQ(((ub - lbv.array()).log() + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
  lp = 0.0;
  EXPECT_MATRIX_EQ(res, stan::math::lub_constrain(input, lb, ub, lp));
  EXPECT_FLOAT_EQ((std::log(ub - lb) + input.array()
                   - 2 * input.array().exp().log1p().array())
                      .sum(),
                  lp);
}

TEST(prob_transform, lub_j) {
  double lp = -17.0;
  double L = 2.0;
  double U = 5.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(L + (U - L) * stan::math::inv_logit(x),
                  stan::math::lub_constrain(x, L, U, lp));
  EXPECT_FLOAT_EQ(-17.0 + log(U - L) + log(stan::math::inv_logit(x))
                      + log(1.0 - stan::math::inv_logit(x)),
                  lp);

  double lp1 = -12.9;
  EXPECT_FLOAT_EQ(-12.9, lp1);

  double lp2 = -19.8;
  double lp2_expected = -19.8;
  EXPECT_FLOAT_EQ(lp2_expected, lp2);

  double lp3 = -422;
  double lp3_expected = -422;
  EXPECT_FLOAT_EQ(lp3_expected, lp3);
}
TEST(prob_transform, lubException) {
  using stan::math::lub_constrain;
  EXPECT_THROW(lub_constrain(5.0, 1.0, 1.0), std::domain_error);
  EXPECT_NO_THROW(lub_constrain(5.0, 1.0, 1.01));
  double lp = 12;
  EXPECT_THROW(lub_constrain(5.0, 1.0, 1.0, lp), std::domain_error);
  EXPECT_NO_THROW(lub_constrain(5.0, 1.0, 1.01, lp));
}
TEST(prob_transform, lub_f) {
  double L = -10.0;
  double U = 27.0;
  double y = 3.0;
  EXPECT_FLOAT_EQ(stan::math::logit((y - L) / (U - L)),
                  stan::math::lub_free(y, L, U));
}
TEST(prob_transform, lub_f_exception) {
  double L = -10.0;
  double U = 27.0;
  EXPECT_THROW(stan::math::lub_free(L - 0.01, L, U), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(U + 0.01, L, U), std::domain_error);

  EXPECT_THROW(stan::math::lub_free((L + U) / 2, U, L), std::domain_error);
}
TEST(prob_transform, lub_rt) {
  double x = -1.0;
  double xc = stan::math::lub_constrain(x, 2.0, 4.0);
  double xcf = stan::math::lub_free(xc, 2.0, 4.0);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::lub_constrain(xcf, 2.0, 4.0);
  EXPECT_FLOAT_EQ(xc, xcfc);
}

TEST(prob_transform, underflow) {
  EXPECT_EQ(0, stan::math::lub_constrain(-1000, 0, 1));
  double lp = 0;
  EXPECT_EQ(0, stan::math::lub_constrain(-1000, 0, 1, lp));
}

TEST(prob_transform, lub_constrain_matrix) {
  Eigen::VectorXd x(4);
  x << -1.0, 1.1, 3.0, 5.0;
  Eigen::VectorXd ub(4);
  ub << stan::math::INFTY, stan::math::INFTY, 6.0, 7.0;
  Eigen::VectorXd lb(4);
  lb << 2.0, stan::math::NEGATIVE_INFTY, stan::math::NEGATIVE_INFTY, 2.0;

  double ubd = 8.0;
  double lbd = -2.0;

  Eigen::VectorXd ub_bad(3);
  Eigen::VectorXd lb_bad(3);

  // matrix, real, real
  {
    Eigen::VectorXd result = stan::math::lub_constrain(x, lbd, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lub_constrain(x(i), lbd, ubd));
    }
    auto x_free = stan::math::lub_free(result, lbd, ubd);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
  }
  // matrix, matrix, real
  {
    Eigen::VectorXd result = stan::math::lub_constrain(x, lb, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lub_constrain(x(i), lb(i), ubd));
    }
    auto x_free = stan::math::lub_free(result, lb, ubd);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lb_bad, ubd),
                 std::invalid_argument);
  }
  // matrix, real, matrix
  {
    Eigen::VectorXd result = stan::math::lub_constrain(x, lbd, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lub_constrain(x(i), lbd, ub(i)));
    }
    auto x_free = stan::math::lub_free(result, lbd, ub);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lbd, ub_bad),
                 std::invalid_argument);
  }
  // matrix, matrix, matrix
  {
    Eigen::VectorXd result = stan::math::lub_constrain(x, lb, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i), stan::math::lub_constrain(x(i), lb(i), ub(i)));
    }
    auto x_free = stan::math::lub_free(result, lb, ub);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lb, ub_bad),
                 std::invalid_argument);
  }
  // matrix, real, real, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::lub_constrain(x, lbd, ubd, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i),
                      stan::math::lub_constrain(x(i), lbd, ubd, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    auto x_free = stan::math::lub_free(result, lbd, ubd);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
  }
  // matrix, matrix, real, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::lub_constrain(x, lb, ubd, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i),
                      stan::math::lub_constrain(x(i), lb(i), ubd, lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    auto x_free = stan::math::lub_free(result, lb, ubd);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lb_bad, ubd),
                 std::invalid_argument);
  }
  // matrix, real, matrix, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::lub_constrain(x, lbd, ub, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i),
                      stan::math::lub_constrain(x(i), lbd, ub(i), lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    auto x_free = stan::math::lub_free(result, lbd, ub);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lbd, ub_bad),
                 std::invalid_argument);
  }
  // matrix, matrix, matrix, lp
  {
    double lp0 = 0.0;
    double lp1 = 0.0;
    Eigen::VectorXd result = stan::math::lub_constrain(x, lb, ub, lp0);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result(i),
                      stan::math::lub_constrain(x(i), lb(i), ub(i), lp1));
    }
    EXPECT_FLOAT_EQ(lp0, lp1);
    auto x_free = stan::math::lub_free(result, lb, ub);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x.coeff(i), x_free.coeff(i));
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lb, ub_bad),
                 std::invalid_argument);
  }
}

TEST(prob_transform, lub_constrain_stdvec) {
  std::vector<double> x{-1.0, 1.1, 3.0, 5.0};
  std::vector<double> ub{stan::math::INFTY, stan::math::INFTY, 6.0, 7.0};
  std::vector<double> lb{2.0, stan::math::NEGATIVE_INFTY,
                         stan::math::NEGATIVE_INFTY, 2.0};

  double ubd = 8.0;
  double lbd = -2.0;
  // These are for dimension checks
  std::vector<double> lb_bad{-2, -3};
  std::vector<double> ub_bad{8, 9};

  // array[], real, real
  {
    auto result = stan::math::lub_constrain(x, lbd, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::lub_constrain(x[i], lbd, ubd));
    }
    auto x_free = stan::math::lub_free(result, lbd, ubd);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
  }
  // array[], array[], real
  {
    auto result = stan::math::lub_constrain(x, lb, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::lub_constrain(x[i], lb[i], ubd));
    }
    auto x_free = stan::math::lub_free(result, lb, ubd);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lb_bad, ubd),
                 std::invalid_argument);
  }
  // array[], real, array[]
  {
    auto result = stan::math::lub_constrain(x, lbd, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::lub_constrain(x[i], lbd, ub[i]));
    }
    auto x_free = stan::math::lub_free(result, lbd, ub);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lbd, ub_bad),
                 std::invalid_argument);
  }
  // array[], array[], array[]
  {
    auto result = stan::math::lub_constrain(x, lb, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      EXPECT_FLOAT_EQ(result[i], stan::math::lub_constrain(x[i], lb[i], ub[i]));
    }
    auto x_free = stan::math::lub_free(result, lb, ub);
    for (size_t i = 0; i < x.size(); ++i) {
      EXPECT_FLOAT_EQ(x[i], x_free[i]);
    }
    EXPECT_THROW(stan::math::lub_constrain(x, lb, ub_bad),
                 std::invalid_argument);
  }
}

TEST(prob_transform, lub_constrain_stdvec_matrix) {
  Eigen::VectorXd x(4);
  x << -1.0, 1.1, 3.0, 5.0;
  Eigen::VectorXd ub(4);
  ub << stan::math::INFTY, stan::math::INFTY, 6.0, 7.0;
  Eigen::VectorXd lb(4);
  lb << 2.0, stan::math::NEGATIVE_INFTY, stan::math::NEGATIVE_INFTY, 2.0;

  double ubd = 8.0;
  double lbd = -2.0;

  Eigen::VectorXd ub_bad(3);
  Eigen::VectorXd lb_bad(3);

  std::vector<Eigen::VectorXd> x_vec{x, x};
  std::vector<Eigen::VectorXd> ub_vec{ub, ub};
  std::vector<Eigen::VectorXd> lb_vec{lb, lb};
  std::vector<Eigen::VectorXd> ub_bad_vec{ub_bad, ub_bad, ub_bad};
  std::vector<Eigen::VectorXd> lb_bad_vec{lb_bad, lb_bad, ub_bad};

  // array[] matrix, array[] matrix, array[] matrix
  {
    auto result = stan::math::lub_constrain(x_vec, lb_vec, ub_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(
            result[i](j),
            stan::math::lub_constrain(x_vec[i](j), lb_vec[i](j), ub_vec[i](j)));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lb_vec, ub_vec);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_bad_vec, ub_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_bad_vec, ub_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_vec, ub_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_bad_vec, ub_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_bad_vec, ub_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_vec, ub_bad_vec),
                 std::invalid_argument);
  }
  // array[] matrix, array[] matrix, real
  {
    auto result = stan::math::lub_constrain(x_vec, lb_vec, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::lub_constrain(
                                          x_vec[i](j), lb_vec[i](j), ubd));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lb_vec, ubd);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_bad_vec, ubd),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_bad_vec, ubd),
                 std::invalid_argument);
  }
  // array[] matrix, real, array[] matrix
  {
    auto result = stan::math::lub_constrain(x_vec, lbd, ub_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::lub_constrain(
                                          x_vec[i](j), lbd, ub_vec[i](j)));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lbd, ub_vec);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lbd, ub_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lbd, ub_bad_vec),
                 std::invalid_argument);
  }
  // array[] matrix, array[] matrix, matrix
  {
    auto result = stan::math::lub_constrain(x_vec, lb_vec, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::lub_constrain(
                                          x_vec[i](j), lb_vec[i](j), ub(j)));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lb_vec, ub);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_vec, ub_bad),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_vec, ub_bad),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_bad_vec, ub),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_bad_vec, ub),
                 std::invalid_argument);
  }
  // array[] matrix, matrix, array[] matrix
  {
    auto result = stan::math::lub_constrain(x_vec, lb, ub_vec);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j), stan::math::lub_constrain(
                                          x_vec[i](j), lb(j), ub_vec[i](j)));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lb, ub_vec);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_bad, ub_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_bad, ub_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb, ub_bad_vec),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb, ub_bad_vec),
                 std::invalid_argument);
  }
  // array[] matrix, matrix, matrix
  {
    auto result = stan::math::lub_constrain(x_vec, lb, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::lub_constrain(x_vec[i](j), lb(j), ub(j)));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lb, ub);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb_bad, ub),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb_bad, ub),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_constrain(x_vec, lb, ub_bad),
                 std::invalid_argument);
    EXPECT_THROW(stan::math::lub_free(x_vec, lb, ub_bad),
                 std::invalid_argument);
  }
  // The rest of the real versions of these don't need the throw test
  // array[] matrix, matrix, real
  {
    auto result = stan::math::lub_constrain(x_vec, lb, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::lub_constrain(x_vec[i](j), lb(j), ubd));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lb, ubd);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
  }
  // array[] matrix, real, matrix
  {
    auto result = stan::math::lub_constrain(x_vec, lbd, ub);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::lub_constrain(x_vec[i](j), lbd, ub(j)));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lbd, ub);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
  }
  // array[] matrix, real, real
  {
    auto result = stan::math::lub_constrain(x_vec, lbd, ubd);
    for (size_t i = 0; i < result.size(); ++i) {
      for (size_t j = 0; j < result[i].size(); ++j) {
        EXPECT_FLOAT_EQ(result[i](j),
                        stan::math::lub_constrain(x_vec[i](j), lbd, ubd));
      }
    }
    auto x_vec_free = stan::math::lub_free(result, lbd, ubd);
    for (size_t i = 0; i < x_vec.size(); ++i) {
      for (size_t j = 0; j < x_vec[i].size(); ++j) {
        EXPECT_FLOAT_EQ(x_vec[i].coeff(j), x_vec_free[i].coeff(j));
      }
    }
  }
}

TEST(prob_transform, lub_free_vec_exception) {
  double ubd = 2.0;
  Eigen::VectorXd ub(2);
  ub << 1.0, 2.0;
  Eigen::VectorXd ub2(2);
  ub2 << -1.0, 5.0;
  std::vector<Eigen::VectorXd> ub_vec = {ub, ub2};

  double lbd = 1.0;
  Eigen::VectorXd lb(2);
  lb << 0.0, 1.0;
  Eigen::VectorXd lb2(2);
  lb2 << -2.0, 4.0;
  std::vector<Eigen::VectorXd> lb_vec = {lb, lb2};

  double xd_bad = 5.0;
  Eigen::VectorXd x_bad(2);
  x_bad << 10.0, 5.0;
  Eigen::VectorXd x_bad2(2);
  x_bad2 << 10.0, 5.0;
  std::vector<Eigen::VectorXd> x_bad_vec = {x_bad, x_bad2};

  EXPECT_THROW(stan::math::lub_free(xd_bad, lbd, ubd), std::domain_error);

  EXPECT_THROW(stan::math::lub_free(x_bad, lbd, ubd), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad, lbd, ub), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad, lb, ubd), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad, lb, ub), std::domain_error);

  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lbd, ubd), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lbd, ub), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lbd, ub_vec), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lb, ubd), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lb, ub), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lb, ub_vec), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lb_vec, ubd), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lb_vec, ub), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(x_bad_vec, lb_vec, ub_vec),
               std::domain_error);
}
