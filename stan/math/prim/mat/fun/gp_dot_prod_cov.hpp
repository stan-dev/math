#ifndef STAN_MATH_PRIM_MAT_FUN_COV_DOT_PROD_HPP
#define STAN_MATH_PRIM_MAT_FUN_COV_DOT_PROD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat/fun/dot_self.hpp>
#include <stan/math/prim/mat/fun/multiply.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns a dot product covariance matrix. A member of Stan's Gaussian
 * Process Library.
 *
 * \f$k(x,x') = \sigma^2 x \cdot x'\f$
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma_squared scalar type of sigma_squared
 *
 * @param x std::vector of elements that can be used in transpose
 *   and multiply
 *    This function assumes each element of x is the same size.
 * @param sigma_squared variance
 * @return dot product covariance matrix that is positive semi-definite
 * @throw std::domain_error if sigma < 0, nan, inf or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma_squared>
Eigen::Matrix<return_type_t<T_x, T_sigma_squared>, Eigen::Dynamic,
              Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<T_x> &x,
                const T_sigma_squared &sigma_squared) {
  check_not_nan("gp_dot_prod_cov", "sigma_squared", sigma_squared);
  check_positive_finite("gp_dot_prod_cov", "sigma_squared", sigma_squared);

  size_t x_size = x.size();
  check_not_nan("gp_dot_prod_cov", "x", x);
  check_finite("gp_dot_prod_cov", "x", x);

  Eigen::Matrix<return_type_t<T_x, T_sigma_squared>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }

  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = sigma_squared * x[i] * x[i];
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = sigma_squared * x[i] * x[j];
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = sigma_squared * x[x_size - 1] * x[x_size - 1];
  return cov;
}

/**
 * Returns a dot product covariance matrix. A member of Stan's Gaussian
 * Process Library.
 *
 * \f$k(x,x') = \sigma^2 x \cdot x'\f$
 *
 * @tparam T_x1 type of first std::vector of double
 * @tparam T_x2 type of second std::vector of double
 * @tparam T_sigma_squared scalar type of sigma_squared
 *
 * @param x1 std::vector of elements that can be used in dot_product
 * @param x2 std::vector of elements that can be used in dot_product
 * @param sigma_squared variance
 * @return dot product covariance matrix
 * @throw std::domain_error if sigma < 0, nan or inf
 *   or if x1 or x2 are nan or inf
 */
template <typename T_x1, typename T_x2, typename T_sigma_squared>
Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma_squared>, Eigen::Dynamic,
              Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<T_x1> &x1, const std::vector<T_x2> &x2,
                const T_sigma_squared &sigma_squared) {
  check_not_nan("gp_dot_prod_cov", "sigma_squared", sigma_squared);
  check_positive_finite("gp_dot_prod_cov", "sigma_squared", sigma_squared);

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  check_not_nan("gp_dot_prod_cov", "x1", x1);
  check_finite("gp_dot_prod_cov", "x1", x1);
  check_not_nan("gp_dot_prod_cov", "x2", x2);
  check_finite("gp_dot_prod_cov", "x2", x2);

  Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma_squared>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = sigma_squared * x1[i] * x2[j];
    }
  }
  return cov;
}

/**
 * Returns a dot product covariance matrix. A member of Stan's Gaussian Process
 * Library.
 *
 * \f$k(x,x') = x^T sigma\_squared I x'\f$
 *
 * @tparam T_x type of first std::vector of double
 * @tparam T_sigma_squared scalar type of sigma_squared
 *
 * @param x std::vector of elements that can be used in dot_product
 * @param sigma_squared variance
 * @return dot product covariance matrix that is positive semi-definite
 * @throw std::domain_error if sigma < 0, nan, inf or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma_squared>
Eigen::Matrix<return_type_t<T_x, T_sigma_squared>, Eigen::Dynamic,
              Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1>> &x,
                const T_sigma_squared &sigma_squared) {
  check_not_nan("gp_dot_prod_cov", "sigma_squared", sigma_squared);
  check_positive_finite("gp_dot_prod_cov", "sigma_squared", sigma_squared);

  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x", x[i]);
    check_finite("gp_dot_prod_cov", "x", x[i]);
  }

  Eigen::Matrix<return_type_t<T_x, T_sigma_squared>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }

  std::vector<Eigen::Matrix<return_type_t<T_x, T_sigma_squared>, Eigen::Dynamic, 1>> x_new(x_size);

  for (size_t i = 0; i < x_size; ++i) {
    x_new[i] = sigma_squared * x[i];
  }

  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = dot_product(x[i], x_new[i]);
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = dot_product(x[i], x_new[j]);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = dot_product(x[x_size - 1], x_new[x_size - 1]);
  return cov;
}

/**
 * Returns a dot product covariance matrix of differing
 * x's. A member of Stan's Gaussian Process Library.
 *
 * \f$k(x,x') = x^T sigma\_squared I x'\f$
 *
 * @tparam T_x1 type of first std::vector of double
 * @tparam T_x2 type of second std::vector of double
 * @tparam T_sigma_squared scalar type of sigma_squared
 *
 * @param x1 std::vector of elements that can be used in dot_product
 * @param x2 std::vector of elements that can be used in dot_product
 * @param sigma_squared variance
 * @return dot product covariance matrix
 * @throw std::domain_error if sigma < 0, nan or inf
 *   or if x1 or x2 are nan or inf
 */
template <typename T_x1, typename T_x2, typename T_sigma_squared>
Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma_squared>, Eigen::Dynamic,
              Eigen::Dynamic>
gp_dot_prod_cov(const std::vector<Eigen::Matrix<T_x1, Eigen::Dynamic, 1>> &x1,
                const std::vector<Eigen::Matrix<T_x2, Eigen::Dynamic, 1>> &x2,
                const T_sigma_squared &sigma_squared) {
  check_not_nan("gp_dot_prod_cov", "sigma_squared", sigma_squared);
  check_positive_finite("gp_dot_prod_cov", "sigma_squared", sigma_squared);

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  for (size_t i = 0; i < x1_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x1", x1[i]);
    check_finite("gp_dot_prod_cov", "x1", x1[i]);
  }
  for (size_t i = 0; i < x2_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x2", x2[i]);
    check_finite("gp_dot_prod_cov", "x2", x2[i]);
  }
  Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma_squared>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  check_size_match("gp_dot_prod_cov", "dimension of elements of x1",
                   size_of(x1[0]), "dimension of elements of x2",
                   size_of(x2[0]));

  std::vector<Eigen::Matrix<return_type_t<T_x2, T_sigma_squared>, Eigen::Dynamic, 1>> x2_new(x2_size);

  for (size_t i = 0; i < x2_size; ++i) {
    x2_new[i] = sigma_squared * x2[i];
  }

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = dot_product(x1[i], x2_new[j]);
    }
  }
  return cov;
}

/**
 * Returns a dot product covariance matrix. A member of Stan's Gaussian Process
 * Library.
 *
 * \f$k(x,x') = x^T diagonal\_matrix(diagonal\_Sigma) x'\f$
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma scalar type of diagonal_Sigma elements
 *
 * @param x std::vector of elements that can be used in dot product.
 *    This function assumes each element of x is the same size.
 * @param diagonal_Sigma vector that is interpreted as a diagonal
 *   covariance matrix
 * @return dot product covariance matrix that is positive semi-definite
 * @throw std::domain_error if sigma < 0, nan, inf or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma>
Eigen::Matrix<return_type_t<T_x, T_sigma>, Eigen::Dynamic, Eigen::Dynamic>
gp_dot_prod_cov(
    const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1>> &x,
    const Eigen::Matrix<T_sigma, Eigen::Dynamic, 1> &diagonal_Sigma) {
  check_not_nan("gp_dot_prod_cov", "diagonal_Sigma", diagonal_Sigma);
  check_positive_finite("gp_dot_prod_cov", "diagonal_Sigma", diagonal_Sigma);

  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x", x[i]);
    check_finite("gp_dot_prod_cov", "x", x[i]);
  }

  Eigen::Matrix<return_type_t<T_x, T_sigma>, Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }

  check_size_match("gp_dot_prod_cov", "dimension of elements of x",
                   size_of(x[0]), "dimension of diagonal_Sigma",
                   diagonal_Sigma.rows());

  std::vector<Eigen::Matrix<return_type_t<T_x, T_sigma>, Eigen::Dynamic, 1>> x_new(x_size);

  for (size_t i = 0; i < x_size; ++i) {
    x_new[i].resize(x[i].rows(), 1);
    for (size_t j = 0; j < diagonal_Sigma.rows(); ++j) {
      x_new[i][j] = diagonal_Sigma(j) * x[i][j];
    }
  }

  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = dot_product(x[i], x_new[i]);
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = dot_product(x[i], x_new[j]);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = dot_product(x[x_size - 1], x_new[x_size - 1]);
  return cov;
}

/**
 * Returns a dot product covariance matrix of differing
 * x's. A member of Stan's Gaussian Process Library.
 *
 * \f$k(x,x') = x^T diagonal\_matrix(diagonal\_Sigma) x'\f$
 *
 * @tparam T_x1 scalar type of x1
 * @tparam T_x2 scalar type of x2
 * @tparam T_sigma scalar type of diagonal_Sigma elements
 *
 * @param x1 std::vector of elements that can be used in dot_product
 * @param x2 std::vector of elements that can be used in dot_product
 * @param diagonal_Sigma vector that is interpreted as a diagonal
 *   covariance matrix
 * @return dot product covariance matrix
 * @throw std::domain_error if sigma < 0, nan or inf
 *   or if x1 or x2 are nan or inf
 */
template <typename T_x1, typename T_x2, typename T_sigma>
Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma>, Eigen::Dynamic,
              Eigen::Dynamic>
gp_dot_prod_cov(
    const std::vector<Eigen::Matrix<T_x1, Eigen::Dynamic, 1>> &x1,
    const std::vector<Eigen::Matrix<T_x2, Eigen::Dynamic, 1>> &x2,
    const Eigen::Matrix<T_sigma, Eigen::Dynamic, 1> &diagonal_Sigma) {
  check_not_nan("gp_dot_prod_cov", "diagonal_Sigma", diagonal_Sigma);
  check_positive_finite("gp_dot_prod_cov", "diagonal_Sigma", diagonal_Sigma);

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  for (size_t i = 0; i < x1_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x1", x1[i]);
    check_finite("gp_dot_prod_cov", "x1", x1[i]);
  }
  for (size_t i = 0; i < x2_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x2", x2[i]);
    check_finite("gp_dot_prod_cov", "x2", x2[i]);
  }
  Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  check_size_match("gp_dot_prod_cov", "dimension of elements of x1",
                   size_of(x1[0]), "dimension of elements of x2",
                   size_of(x2[0]));

  check_size_match("gp_dot_prod_cov", "dimension of elements of x1",
                   size_of(x1[0]), "dimension of Sigma", diagonal_Sigma.rows());

  check_size_match("gp_dot_prod_cov", "dimension of elements of x2",
                   size_of(x2[0]), "dimension of Sigma", diagonal_Sigma.rows());

  std::vector<Eigen::Matrix<return_type_t<T_x2, T_sigma>, Eigen::Dynamic, 1>> x2_new(x2_size);

  for (size_t i = 0; i < x2_size; ++i) {
    x2_new[i].resize(x2[i].rows(), 1);
    for (size_t j = 0; j < diagonal_Sigma.rows(); ++j) {
      x2_new[i][j] = diagonal_Sigma[j] * x2[i][j];
    }
  }

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = dot_product(x1[i], x2_new[j]);
    }
  }
  return cov;
}

/**
 * Returns a dot product covariance matrix. A member of Stan's Gaussian Process
 * Library.
 *
 * \f$k(x,x') = x^T \Sigma x'\f$
 *
 * @tparam T_x type of std::vector of elements
 * @tparam T_sigma scalar type of sigma elements
 *
 * @param x std::vector of elements that can be used in dot product.
 *    This function assumes each element of x is the same size.
 * @param Sigma covariance matrix
 * @return dot product covariance matrix that is positive semi-definite
 * @throw std::domain_error if sigma < 0, nan, inf or
 *   x is nan or infinite
 */
template <typename T_x, typename T_sigma>
Eigen::Matrix<return_type_t<T_x, T_sigma>, Eigen::Dynamic, Eigen::Dynamic>
gp_dot_prod_cov(
    const std::vector<Eigen::Matrix<T_x, Eigen::Dynamic, 1>> &x,
    const Eigen::Matrix<T_sigma, Eigen::Dynamic, Eigen::Dynamic> &Sigma) {
  check_finite("gp_dot_prod_cov", "Sigma", Sigma);
  check_pos_definite("gp_dot_prod_cov", "Sigma", Sigma);

  size_t x_size = x.size();
  for (size_t i = 0; i < x_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x", x[i]);
    check_finite("gp_dot_prod_cov", "x", x[i]);
  }

  Eigen::Matrix<return_type_t<T_x, T_sigma>, Eigen::Dynamic, Eigen::Dynamic>
      cov(x_size, x_size);
  if (x_size == 0) {
    return cov;
  }

  check_size_match("gp_dot_prod_cov", "dimension of elements of x",
                   size_of(x[0]), "dimension of Sigma", Sigma.cols());

  std::vector<Eigen::Matrix<return_type_t<T_x, T_sigma>, Eigen::Dynamic, 1>> x_new(x_size);

  for (size_t i = 0; i < x_size; ++i) {
    x_new[i] = multiply(Sigma, x[i]);
  }

  for (size_t i = 0; i < (x_size - 1); ++i) {
    cov(i, i) = dot_product(x[i], x_new[i]);
    for (size_t j = i + 1; j < x_size; ++j) {
      cov(i, j) = dot_product(x[i], x_new[j]);
      cov(j, i) = cov(i, j);
    }
  }
  cov(x_size - 1, x_size - 1) = dot_product(x[x_size - 1], x_new[x_size - 1]);
  return cov;
}

/**
 * Returns a dot product covariance matrix of differing
 * x's. A member of Stan's Gaussian Process Library.
 *
 * \f$k(x,x') = x^T \Sigma x'\f$
 *
 * @tparam T_x1 scalar type of x1
 * @tparam T_x2 scalar type of x2
 * @tparam T_sigma scalar type of Sigma
 *
 * @param x1 std::vector of elements that can be used in dot_product
 * @param x2 std::vector of elements that can be used in dot_product
 * @param Sigma covariance matrix
 * @return dot product covariance matrix
 * @throw std::domain_error if sigma < 0, nan or inf
 *   or if x1 or x2 are nan or inf
 */
template <typename T_x1, typename T_x2, typename T_sigma>
Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma>, Eigen::Dynamic,
              Eigen::Dynamic>
gp_dot_prod_cov(
    const std::vector<Eigen::Matrix<T_x1, Eigen::Dynamic, 1>> &x1,
    const std::vector<Eigen::Matrix<T_x2, Eigen::Dynamic, 1>> &x2,
    const Eigen::Matrix<T_sigma, Eigen::Dynamic, Eigen::Dynamic> &Sigma) {
  check_finite("gp_dot_prod_cov", "Sigma", Sigma);
  check_pos_definite("gp_dot_prod_cov", "Sigma", Sigma);

  size_t x1_size = x1.size();
  size_t x2_size = x2.size();
  for (size_t i = 0; i < x1_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x1", x1[i]);
    check_finite("gp_dot_prod_cov", "x1", x1[i]);
  }
  for (size_t i = 0; i < x2_size; ++i) {
    check_not_nan("gp_dot_prod_cov", "x2", x2[i]);
    check_finite("gp_dot_prod_cov", "x2", x2[i]);
  }
  Eigen::Matrix<return_type_t<T_x1, T_x2, T_sigma>, Eigen::Dynamic,
                Eigen::Dynamic>
      cov(x1_size, x2_size);

  if (x1_size == 0 || x2_size == 0) {
    return cov;
  }

  check_size_match("gp_dot_prod_cov", "dimension of elements of x1",
                   size_of(x1[0]), "dimension of elements of x2",
                   size_of(x2[0]));

  check_size_match("gp_dot_prod_cov", "dimension of elements of x1",
                   size_of(x1[0]), "dimension of Sigma", Sigma.cols());

  check_size_match("gp_dot_prod_cov", "dimension of elements of x2",
                   size_of(x2[0]), "dimension of Sigma", Sigma.cols());

  std::vector<Eigen::Matrix<return_type_t<T_x2, T_sigma>, Eigen::Dynamic, 1>> x2_new(x2_size);

  for (size_t i = 0; i < x2_size; ++i) {
    x2_new[i] = multiply(Sigma, x2[i]);
  }

  for (size_t i = 0; i < x1_size; ++i) {
    for (size_t j = 0; j < x2_size; ++j) {
      cov(i, j) = dot_product(x1[i], x2_new[j]);
    }
  }
  return cov;
}
}  // namespace math
}  // namespace stan
#endif
