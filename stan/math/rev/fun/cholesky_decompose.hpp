#ifndef STAN_MATH_REV_FUN_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_REV_FUN_CHOLESKY_DECOMPOSE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/opencl.hpp>
#endif

#include <algorithm>
#include <vector>

namespace stan {
namespace math {


#ifdef STAN_OPENCL
class cholesky_opencl : public vari {
 public:
  int M_;
  vari** vari_ref_A_;
  vari** vari_ref_L_;

  /**
   * Constructor for OpenCL Cholesky function.
   *
   * Stores varis for A.  Instantiates and stores varis for L.
   * Instantiates and stores dummy vari for upper triangular part of var
   * result returned in cholesky_decompose function call
   *
   * variRefL aren't on the chainable autodiff stack, only used for storage
   * and computation. Note that varis for L are constructed externally in
   * cholesky_decompose.
   *
   * @param A matrix
   * @param L_A Cholesky factor of A
   */
  cholesky_opencl(const Eigen::Matrix<var, -1, -1>& A,
                  const Eigen::Matrix<double, -1, -1>& L_A)
      : vari(0.0),
        M_(A.rows()),
        vari_ref_A_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A.rows() * (A.rows() + 1) / 2)),
        vari_ref_L_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            A.rows() * (A.rows() + 1) / 2)) {
    size_t pos = 0;
    for (size_type j = 0; j < M_; ++j) {
      for (size_type i = j; i < M_; ++i) {
        vari_ref_A_[pos] = A.coeffRef(i, j).vi_;
        vari_ref_L_[pos] = new vari(L_A.coeffRef(i, j), false);
        ++pos;
      }
    }
  }

  /**
   * Symbolic adjoint calculation for Cholesky factor A
   *
   * @param L_val value of Cholesky factor
   * @param L_adj adjoint of Cholesky factor
   */
  inline void symbolic_rev(matrix_cl<double>& L_val, matrix_cl<double>& L_adj) {
    L_adj = transpose(L_val) * L_adj;
    L_adj.triangular_transpose<TriangularMapCL::LowerToUpper>();
    L_val = transpose(tri_inverse(L_val));
    L_adj = L_val * transpose(L_val * L_adj);
    L_adj.triangular_transpose<TriangularMapCL::LowerToUpper>();
  }

  /**
   * Reverse mode differentiation algorithm using OpenCL
   *
   * Reference:
   *
   * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
   *
   */
  virtual void chain() {
    const int packed_size = M_ * (M_ + 1) / 2;
    Eigen::Map<Eigen::Matrix<vari*, Eigen::Dynamic, 1>> L_cpu(
        vari_ref_L_, M_ * (M_ + 1) / 2);
    Eigen::VectorXd L_val_cpu = L_cpu.val();
    Eigen::VectorXd L_adj_cpu = L_cpu.adj();
    matrix_cl<double> L_val = packed_copy<matrix_cl_view::Lower>(L_val_cpu, M_);
    matrix_cl<double> L_adj = packed_copy<matrix_cl_view::Lower>(L_adj_cpu, M_);
    int block_size
        = M_ / opencl_context.tuning_opts().cholesky_rev_block_partition;
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

      auto&& R_val = block(L_val, j, 0, k_j_ind, j);
      auto&& R_adj = block(L_adj, j, 0, k_j_ind, j);
      matrix_cl<double> D_val = block(L_val, j, j, k_j_ind, k_j_ind);
      matrix_cl<double> D_adj = block(L_adj, j, j, k_j_ind, k_j_ind);
      auto&& B_val = block(L_val, k, 0, m_k_ind, j);
      auto&& B_adj = block(L_adj, k, 0, m_k_ind, j);
      auto&& C_val = block(L_val, k, j, m_k_ind, k_j_ind);
      auto&& C_adj = block(L_adj, k, j, m_k_ind, k_j_ind);

      C_adj = C_adj * tri_inverse(D_val);
      B_adj = B_adj - C_adj * R_val;
      D_adj = D_adj - transpose(C_adj) * C_val;

      symbolic_rev(D_val, D_adj);

      R_adj = R_adj - transpose(C_adj) * B_val - D_adj * R_val;
      D_adj = diagonal_multiply(D_adj, 0.5);

      block(L_adj, j, j, k_j_ind, k_j_ind) = D_adj;
    }
    L_adj.view(matrix_cl_view::Lower);
    std::vector<double> L_adj_cpu_res = packed_copy(L_adj);
    for (size_type j = 0; j < packed_size; ++j) {
      vari_ref_A_[j]->adj_ += L_adj_cpu_res[j];
    }
  }
};
#endif


namespace internal {
/**
 * Reverse mode differentiation algorithm reference:
 *
 * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
 *
 */
template <typename T1, typename T2, typename T3>
auto cholesky_lambda(T1 L_A, T2 L, T3 A_ref) {
  return [=]() mutable {
    using Eigen::Lower;
    using Eigen::StrictlyUpper;
    using Eigen::Upper;
    auto L_adj = L.adj().eval();
    const int M_ = L_A.rows();
    int block_size_ = std::max(M_ / 8, 8);
    block_size_ = std::min(block_size_, 128);
    size_t pos = 0;
    for (int k = M_; k > 0; k -= block_size_) {
      int j = std::max(0, k - block_size_);
      auto R = L_A.block(j, 0, k - j, j);
      // This is on purpose, all the other ones are views and we write to this
      Eigen::MatrixXd D = L_A.block(j, j, k - j, k - j);
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
    A_ref.adj().template triangularView<Eigen::Lower>() = L_adj;
  };
}
}

/**
 * Reverse mode specialization of Cholesky decomposition
 *
 * Internally calls Eigen::LLT rather than using
 * stan::math::cholesky_decompose in order to use an inplace decomposition.
 *
 * Note chainable stack varis are created below in Matrix<var, -1, -1>
 *
 * @param A Matrix
 * @return L Cholesky factor of A
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline Eigen::Matrix<var, T::RowsAtCompileTime, T::ColsAtCompileTime>
cholesky_decompose(const T& A) {
  arena_t<T> A_ref = A;
  arena_t<promote_scalar_t<double, T>> L_A(value_of(A_ref));
  check_not_nan("cholesky_decompose", "A", L_A);
#ifdef STAN_OPENCL
  L_A = cholesky_decompose(L_A);
#else
  check_symmetric("cholesky_decompose", "A", L_A);
  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>, Eigen::Lower> L_factor(L_A);
  check_pos_definite("cholesky_decompose", "A", L_factor);
#endif
  // Memory allocated in arena.
  // cholesky_scalar gradient faster for small matrices compared to
  // cholesky_block
  vari* dummy = new vari(0.0, false);
  arena_t<T> L(L_A.rows(), L_A.cols());
  for (Eigen::Index j = 0; j < L_A.rows(); ++j) {
    for (Eigen::Index i = 0; i < L_A.rows(); ++i) {
      if (j > i) {
        L.coeffRef(i, j) = dummy;
      } else {
        L.coeffRef(i, j) = new vari(L_A.coeffRef(i, j), false);
      }
    }
  }
  if (L_A.rows() <= 35) {
    reverse_pass_callback([L_A, L, A_ref]() mutable {
      size_t N = A_ref.rows();
      int pos = (N * (N + 1) / 2) - 1;
      auto adjL = L.adj().eval();
      // TODO(Steve): Write this in column major order
      for (int i = N - 1; i >= 0; --i) {
        for (int j = i; j >= 0; --j) {
          if (i == j) {
            A_ref.adj()(i, j) = 0.5 * adjL.coeff(i, j) / L_A.coeff(i, j);
          } else {
            A_ref.adj()(i, j) = adjL.coeff(i, j) / L_A.coeff(j, j);
            adjL.coeffRef(j, j)
                -= adjL.coeff(i, j) * L_A.coeff(i, j) / L_A.coeff(j, j);
          }
          for (int k = j - 1; k >= 0; --k) {
            adjL.coeffRef(i, k) -= A_ref.adj().coeff(i, j) * L_A.coeff(j, k);
            adjL.coeffRef(j, k) -= A_ref.adj().coeff(i, j) * L_A.coeff(i, k);
          }
        }
      }
    });
  } else {
    #ifdef STAN_OPENCL
        if (L_A.rows()
            > opencl_context.tuning_opts().cholesky_size_worth_transfer) {
          cholesky_opencl* baseVari = new cholesky_opencl(A_ref, L_A);
          internal::set_lower_tri_coeff_ref(L, baseVari->vari_ref_L_);
        } else {
          reverse_pass_callback(internal::cholesky_lambda(L_A, L, A_ref));
        }
    #else
    reverse_pass_callback(internal::cholesky_lambda(L_A, L, A_ref));
    #endif
  }

  return L;
}

/**
 * Reverse mode differentiation algorithm reference:
 *
 * Iain Murray: Differentiation of the Cholesky decomposition, 2016.
 *
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto
cholesky_decompose(const T& A) {
  check_symmetric("cholesky_decompose", "A", A.val());
  check_not_nan("cholesky_decompose", "A", A.val());
  arena_t<T> L = A.val().llt().matrixL();
  if (A.rows() <= 35) {
    reverse_pass_callback([L, A]() mutable {
      const size_t N = A.rows();
      auto adjL = L.adj().eval();
      // TODO(Steve): Write this in column major order
      for (int i = N - 1; i >= 0; --i) {
        for (int j = i; j >= 0; --j) {
          if (i == j) {
            A.adj().coeffRef(i, j) = 0.5 * adjL.coeff(i, j) / L.val().coeff(i, j);
          } else {
            A.adj().coeffRef(i, j) = adjL.coeff(i, j) / L.val().coeffRef(j, j);
            adjL.coeffRef(j, j)
                -= adjL.coeffRef(i, j) * L.val().coeffRef(i, j) / L.val().coeffRef(j, j);
          }
          for (int k = j - 1; k >= 0; --k) {
            adjL.coeffRef(i, k) -= A.adj().coeff(i, j) * L.val().coeffRef(j, k);
            adjL.coeffRef(j, k) -= A.adj().coeff(i, j) * L.val().coeffRef(i, k);
          }
        }
      }
    });
  } else {
    reverse_pass_callback([L, A]() mutable {
      using Eigen::Block;
      using Eigen::Lower;
      using Eigen::MatrixXd;
      using Eigen::StrictlyUpper;
      using Eigen::Upper;
      using Block_ = Eigen::Block<Eigen::MatrixXd>;

      const int M_ = L.val().rows();
      auto L_adj = L.adj().eval();
      int block_size_ = std::max(M_ / 8, 8);
      block_size_ = std::min(block_size_, 128);
      for (int k = M_; k > 0; k -= block_size_) {
        int j = std::max(0, k - block_size_);
        auto R = L.val().block(j, 0, k - j, j);
        // This is on purpose, all the other ones are views and we write to this
        Eigen::MatrixXd D = L.val().block(j, j, k - j, k - j);
        auto B = L.val().block(k, 0, M_ - k, j);
        auto C = L.val().block(k, j, M_ - k, k - j);
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
      A.adj().template triangularView<Eigen::Lower>() = L.adj();
    });
  }
  return L;
}

}  // namespace math
}  // namespace stan
#endif
