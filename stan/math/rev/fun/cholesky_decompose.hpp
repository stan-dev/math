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

#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#endif

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
    adjL.template triangularView<Eigen::Lower>() = L.adj();
    for (int i = N - 1; i >= 0; --i) {
      for (int j = i; j >= 0; --j) {
        if (i == j) {
          A.adj()(i, j) = 0.5 * adjL.coeff(i, j) / L_A.coeff(i, j);
        } else {
          A.adj()(i, j) = adjL.coeff(i, j) / L_A.coeff(j, j);
          adjL.coeffRef(j, j)
              -= adjL.coeff(i, j) * L_A.coeff(i, j) / L_A.coeff(j, j);
        }
        for (int k = j - 1; k >= 0; --k) {
          adjL.coeffRef(i, k) -= A.adj().coeff(i, j) * L_A.coeff(j, k);
          adjL.coeffRef(j, k) -= A.adj().coeff(i, j) * L_A.coeff(i, k);
        }
      }
    }
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
    using Block_ = Eigen::Block<Eigen::MatrixXd>;
    using Eigen::Lower;
    using Eigen::StrictlyUpper;
    using Eigen::Upper;
    Eigen::MatrixXd L_adj(L.rows(), L.cols());
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
    A.adj().template triangularView<Eigen::Lower>() = L_adj;
  };
}

#ifdef STAN_OPENCL
/**
 * Reverse mode differentiation for Cholesky using OpenCL
 * Reverse mode differentiation algorithm reference:
 *
 * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
 *
 */
template <typename AMat, typename LVari>
inline auto opencl_cholesky_lambda(AMat& arena_A, LVari& vari_L) {
  return [arena_A, vari_L]() mutable {
    const int M_ = arena_A.rows();
    const int packed_size = M_ * (M_ + 1) / 2;
    Eigen::Map<Eigen::Matrix<vari*, Eigen::Dynamic, 1>> L_cpu(
        vari_L, M_ * (M_ + 1) / 2);
    Eigen::VectorXd L_val_cpu = L_cpu.val();
    Eigen::VectorXd L_adj_cpu = L_cpu.adj();
    matrix_cl<double> L_val = packed_copy<matrix_cl_view::Lower>(L_val_cpu, M_);
    matrix_cl<double> L_adj = packed_copy<matrix_cl_view::Lower>(L_adj_cpu, M_);
    int block_size
        = M_
          / std::max(1,
                     opencl_context.tuning_opts().cholesky_rev_block_partition);
    block_size = std::max(block_size, 8);
    block_size = std::min(
        block_size, opencl_context.tuning_opts().cholesky_rev_min_block_size);
    // The following is an OpenCL implementation of
    // the chain() function from the cholesky_block
    // vari class implementation
    for (int k = M_; k > 0; k -= block_size) {
      const int j = std::max(0, k - block_size);
      const int k_j_ind = k - j;
      const int m_k_ind = M_ - k;

      auto&& R_val = block_zero_based(L_val, j, 0, k_j_ind, j);
      auto&& R_adj = block_zero_based(L_adj, j, 0, k_j_ind, j);
      matrix_cl<double> D_val = block_zero_based(L_val, j, j, k_j_ind, k_j_ind);
      matrix_cl<double> D_adj = block_zero_based(L_adj, j, j, k_j_ind, k_j_ind);
      auto&& B_val = block_zero_based(L_val, k, 0, m_k_ind, j);
      auto&& B_adj = block_zero_based(L_adj, k, 0, m_k_ind, j);
      auto&& C_val = block_zero_based(L_val, k, j, m_k_ind, k_j_ind);
      auto&& C_adj = block_zero_based(L_adj, k, j, m_k_ind, k_j_ind);

      C_adj = C_adj * tri_inverse(D_val);
      B_adj = B_adj - C_adj * R_val;
      D_adj = D_adj - transpose(C_adj) * C_val;

      D_adj = transpose(D_val) * D_adj;
      D_adj.triangular_transpose<TriangularMapCL::LowerToUpper>();
      D_val = transpose(tri_inverse(D_val));
      D_adj = D_val * transpose(D_val * D_adj);
      D_adj.triangular_transpose<TriangularMapCL::LowerToUpper>();

      R_adj = R_adj - transpose(C_adj) * B_val - D_adj * R_val;
      diagonal(D_adj) = diagonal(D_adj) * 0.5;

      block_zero_based(L_adj, j, j, k_j_ind, k_j_ind) = D_adj;
    }
    L_adj.view(matrix_cl_view::Lower);
    std::vector<double> L_adj_cpu_res = packed_copy(L_adj);
    int pos = 0;
    for (Eigen::Index j = 0; j < arena_A.rows(); ++j) {
      for (Eigen::Index i = j; i < arena_A.rows(); ++i) {
        arena_A.coeffRef(i, j)->adj_ += L_adj_cpu_res[pos];
        pos++;
      }
    }
  };
}

#endif
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
#ifdef STAN_OPENCL
  L_A = cholesky_decompose(L_A);
#else
  check_symmetric("cholesky_decompose", "A", A);
  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower> L_factor(L_A);
  check_pos_definite("cholesky_decompose", "m", L_factor);
#endif
  L_A.template triangularView<Eigen::StrictlyUpper>().setZero();
  // looping gradient calcs faster for small matrices compared to
  // cholesky_block
  vari* dummy = new vari(0.0, false);
  arena_t<EigMat> L(L_A.rows(), L_A.cols());
  if (L_A.rows() <= 35) {
    internal::initialize_return(L, L_A, dummy);
    reverse_pass_callback(internal::unblocked_cholesky_lambda(L_A, L, arena_A));
  } else {
#ifdef STAN_OPENCL
    if (L_A.rows()
        > opencl_context.tuning_opts().cholesky_size_worth_transfer) {
      vari** vari_L = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
          arena_A.rows() * (arena_A.rows() + 1) / 2);
      size_t pos = 0;
      for (Eigen::Index j = 0; j < arena_A.rows(); ++j) {
        for (Eigen::Index i = j; i < arena_A.rows(); ++i) {
          vari_L[pos] = new vari(L_A.coeffRef(i, j), false);
          L.coeffRef(i, j).vi_ = vari_L[pos];
          ++pos;
        }
        for (Eigen::Index k = 0; k < j; ++k) {
          L.coeffRef(k, j).vi_ = dummy;
        }
      }
      reverse_pass_callback(internal::opencl_cholesky_lambda(arena_A, vari_L));
    } else {
      internal::initialize_return(L, L_A, dummy);
      reverse_pass_callback(internal::cholesky_lambda(L_A, L, arena_A));
    }
#else
    internal::initialize_return(L, L_A, dummy);
    reverse_pass_callback(internal::cholesky_lambda(L_A, L, arena_A));
#endif
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
