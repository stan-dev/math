#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/err/check_pos_definite.hpp>
#include <stan/math/prim/err/check_square.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <algorithm>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename LMat, typename LAMat>
inline void initialize_return(LMat& L, const LAMat& L_A, vari*& dummy) {
  for (Eigen::Index j = 0; j < L_A.rows(); ++j) {
    for (Eigen::Index i = 0; i < L_A.rows(); ++i) {
      if (j > i) {
        L.coeffRef(i, j) = dummy;
      } else {
        L.coeffRef(i, j) = new vari(L_A.coeffRef(i, j), false);
      }
    }
  }
}

/**
 * Reverse mode differentiation algorithm reference:
 *
 * Mike Giles. An extended collection of matrix derivative results for
 * forward and reverse mode AD.  Jan. 2008.
 *
 * Note algorithm  as laid out in Giles is row-major, so Eigen::Matrices
 * are explicitly storage order RowMajor, whereas Eigen defaults to
 * ColumnMajor. Also note algorithm starts by calculating the adjoint for
 * A(M_ - 1, M_ - 1), hence pos on line 94 is decremented to start at pos
 * = M_ * (M_ + 1) / 2.
 */
template <typename T1, typename T2, typename T3>
inline auto unblocked_cholesky_lambda(T1& L_A, T2& L, T3& A) {
  return [L_A, L, A]() mutable {
    const size_t N = A.rows();
    // Algorithm is in rowmajor so we make the adjoint copy rowmajor
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> adjL(L.rows(), L.cols());
    Eigen::MatrixXd adjA = Eigen::MatrixXd::Zero(L.rows(), L.cols());
    adjL.template triangularView<Eigen::Lower>() = L.adj();
    for (int i = N - 1; i >= 0; --i) {
      for (int j = i; j >= 0; --j) {
        if (i == j) {
          adjA.coeffRef(i, j) = 0.5 * adjL.coeff(i, j) / L_A.coeff(i, j);
        } else {
          adjA.coeffRef(i, j) = adjL.coeff(i, j) / L_A.coeff(j, j);
          adjL.coeffRef(j, j)
              -= adjL.coeff(i, j) * L_A.coeff(i, j) / L_A.coeff(j, j);
        }
        for (int k = j - 1; k >= 0; --k) {
          adjL.coeffRef(i, k) -= adjA.coeff(i, j) * L_A.coeff(j, k);
          adjL.coeffRef(j, k) -= adjA.coeff(i, j) * L_A.coeff(i, k);
        }
      }
    }
    A.adj() += adjA;
  };
}

/**
 * Reverse mode differentiation algorithm reference:
 *
 * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
 *
 */
template <typename T1, typename T2, typename T3>
inline auto cholesky_lambda(T1& L_A, T2& L, T3& A) {
  return [L_A, L, A]() mutable {
    using Eigen::Lower;
    using Eigen::StrictlyUpper;
    using Eigen::Upper;
    Eigen::MatrixXd L_adj = Eigen::MatrixXd::Zero(L.rows(), L.cols());
    L_adj.template triangularView<Eigen::Lower>() = L.adj();
    const int M_ = L_A.rows();
    int block_size_ = std::max(M_ / 8, 8);
    block_size_ = std::min(block_size_, 128);
    for (int k = M_; k > 0; k -= block_size_) {
      int j = std::max(0, k - block_size_);
      auto R = L_A.block(j, 0, k - j, j);
      auto D = L_A.block(j, j, k - j, k - j).eval();
      auto B = L_A.block(k, 0, M_ - k, j);
      auto C = L_A.block(k, j, M_ - k, k - j);
      auto R_adj = L_adj.block(j, 0, k - j, j);
      auto D_adj = L_adj.block(j, j, k - j, k - j);
      auto B_adj = L_adj.block(k, 0, M_ - k, j);
      auto C_adj = L_adj.block(k, j, M_ - k, k - j);
      D.transposeInPlace();
      if (C_adj.size() > 0) {
        C_adj = D.template triangularView<Upper>()
                    .solve(C_adj.transpose())
                    .transpose();
        B_adj.noalias() -= C_adj * R;
        D_adj.noalias() -= C_adj.transpose() * C;
      }
      D_adj = (D * D_adj.template triangularView<Lower>()).eval();
      D_adj.template triangularView<StrictlyUpper>()
          = D_adj.adjoint().template triangularView<StrictlyUpper>();
      D.template triangularView<Upper>().solveInPlace(D_adj);
      D.template triangularView<Upper>().solveInPlace(D_adj.transpose());
      R_adj.noalias() -= C_adj.transpose() * B;
      R_adj.noalias() -= D_adj.template selfadjointView<Lower>() * R;
      D_adj.diagonal() *= 0.5;
    }
    A.adj().template triangularView<Eigen::Lower>() += L_adj;
  };
}
}  // namespace internal

/**
 * Reverse mode specialization of cholesky decomposition
 *
 * Internally calls Eigen::LLT rather than using
 * stan::math::cholesky_decompose in order to use an inplace decomposition.
 *
 * Note chainable stack varis are created below in Matrix<var, -1, -1>
 *
 * @param A Matrix
 * @return L cholesky factor of A
 */
template <typename EigMat, require_eigen_vt<is_var, EigMat>* = nullptr>
inline auto cholesky_decompose(const EigMat& A) {
  check_square("cholesky_decompose", "A", A);
  arena_t<EigMat> arena_A = A;
  arena_t<Eigen::Matrix<double, -1, -1>> L_A(arena_A.val());

  check_symmetric("cholesky_decompose", "A", A);
  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower> L_factor(L_A);
  check_pos_definite("cholesky_decompose", "m", L_factor);

  L_A.template triangularView<Eigen::StrictlyUpper>().setZero();
  // looping gradient calcs faster for small matrices compared to
  // cholesky_block
  vari* dummy = new vari(0.0, false);
  arena_t<EigMat> L(L_A.rows(), L_A.cols());
  if (L_A.rows() <= 35) {
    internal::initialize_return(L, L_A, dummy);
    reverse_pass_callback(internal::unblocked_cholesky_lambda(L_A, L, arena_A));
  } else {
    internal::initialize_return(L, L_A, dummy);
    reverse_pass_callback(internal::cholesky_lambda(L_A, L, arena_A));
  }
  return plain_type_t<EigMat>(L);
}

/**
 * Reverse mode specialization of Cholesky decomposition
 *
 * Internally calls Eigen::LLT rather than using
 * stan::math::cholesky_decompose in order to use an inplace decomposition.
 *
 * Note chainable stack varis are created below in Matrix<var, -1, -1>
 * @tparam T A `var_value` holding an inner eigen type.
 * @param A A square positive definite matrix with no nan values.
 * @return L Cholesky factor of A
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto cholesky_decompose(const T& A) {
  check_symmetric("cholesky_decompose", "A", A.val());
  plain_type_t<T> L = cholesky_decompose(A.val());
  if (A.rows() <= 35) {
    reverse_pass_callback(internal::unblocked_cholesky_lambda(L.val(), L, A));
  } else {
    reverse_pass_callback(internal::cholesky_lambda(L.val(), L, A));
  }
  return L;
}

}  // namespace math
}  // namespace stan
#endif
