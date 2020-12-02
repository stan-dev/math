#ifndef STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_GENERALIZED_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/chol2inv.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/tcrossprod.hpp>
#include <stan/math/prim/fun/crossprod.hpp>
#include <stan/math/prim/fun/quad_form.hpp>
#include <stan/math/prim/fun/inverse_spd.hpp>
#include <stan/math/prim/fun/mdivide_left_spd.hpp>
#include <stan/math/prim/fun/mdivide_right_spd.hpp>
#include <stan/math/prim/fun/inverse.hpp>

namespace stan {
namespace math {

/**
 * Returns the Moore-Penrose generalized inverse of the specified matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam n number of rows, can be Eigen::Dynamic
 * @tparam m number of columns, can be Eigen::Dynamic
 *
 * @param G specified matrix
 * @return Generalized inverse of the matrix (an empty matrix if the specified
 * matrix has size zero).
 */
template <typename EigMat,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
generalized_inverse(const EigMat& G) {
  using value_t = value_type_t<EigMat>;
  if (G.size() == 0) {
    return {};
  }

  const auto n = G.rows();
  const auto m = G.cols();

  if (G.rows() == G.cols()) {
    return inverse(G);
  }

  if (n < m) {
    return transpose(mdivide_left_spd(A, tcrossprod(G)));
  } else {
    return transpose(mdivide_right_spd(G, crossprod(A)));
  }
}

/**
 * Returns the Moore-Penrose generalized inverse of the specified matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam n number of rows, can be Eigen::Dynamic
 * @tparam m number of columns, can be Eigen::Dynamic
 *
 * @param G specified matrix
 * @param a diagonal jitter
 * @return Generalized inverse of the matrix (an empty matrix if the specified
 * matrix has size zero).
 */
template <typename EigMat,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
generalized_inverse(const EigMat& G, double a) {
  using value_t = value_type_t<EigMat>;
  if (G.size() == 0) {
    return {};
  }

  const auto n = G.rows();
  const auto m = G.cols();

  if (G.rows() == G.cols()) {
    return inverse(G);
  }

  if (n < m) {
    Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> A = tcrossprod(G);
    A.diagonal().array()
        += Eigen::Array<double, Eigen::Dynamic, 1>::Constant(n, a);
    return transpose(mdivide_left_spd(A, G));
  } else {
    Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic> A = crossprod(G);
    A.diagonal().array()
        += Eigen::Array<double, Eigen::Dynamic, 1>::Constant(m, a);
    return transpose(mdivide_right_spd(G, A));
  }
}

}  // namespace math
}  // namespace stan

#endif
