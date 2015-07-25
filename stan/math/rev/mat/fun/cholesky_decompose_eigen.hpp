#ifndef STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_EIGEN_HPP
#define STAN_MATH_REV_MAT_FUN_CHOLESKY_DECOMPOSE_EIGEN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_pos_definite.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <numeric>

namespace stan {
  namespace math {

    Eigen::SparseMatrix<double>
    id_p_commutation_matrix(int n) {
      typedef Eigen::Triplet<double> triple;
      std::vector<triple> sparse_list;
      sparse_list.reserve(static_cast<size_t>(n * n));
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
          int rowpos = i * n + j;
          int colpos = j * n + i;
          sparse_list.push_back(triple(rowpos, colpos, 1));
        }
      Eigen::SparseMatrix<double> m(n * n, n * n);
      m.setFromTriplets(sparse_list.begin(), sparse_list.end());
      return m + Eigen::MatrixXd::Identity(n * n, n * n).sparseView();
    }

    template <typename T>
    Eigen::SparseMatrix<T>
    elimination_matrix(int n) {
      typedef Eigen::Triplet<T> triple;
      std::vector<triple> sparse_list;
      int p = n * (n + 1) / 2;
      sparse_list.reserve(static_cast<size_t>(p));
      int groupctr = 0;
      int k = 0;
      int l = 1;
      int nuse = n;
      for (int i = 0; i < p; ++i) {
        ++k;
        if (k > nuse) {
          k = 1;
          --nuse;
          groupctr += l;
          ++l;
        }
        int j = i + groupctr;
        sparse_list.push_back(triple(i, j, 1));
      }
      Eigen::SparseMatrix<T> m(p, n * n);
      m.setFromTriplets(sparse_list.begin(), sparse_list.end());
      return m;
    }

    template <typename T>
    Eigen::SparseMatrix<T>
    chol_kron_eye(Eigen::Matrix<T, -1, -1> L) {
      typedef Eigen::Triplet<T> triple;
      std::vector<triple> sparse_list;
      int n = L.rows();
      int p = n * (n + 1) / 2;
      sparse_list.reserve(static_cast<size_t>(p * n));
      for (int j = 0; j < n; ++j)
        for (int i = j; i < n; ++i) {
          int row = n * i;
          int col = n * j;
          for (int m = 0; m < L.rows(); ++m)
            sparse_list.push_back(triple(row + m, col + m, L(i, j)));
        }
      Eigen::SparseMatrix<T> m(n * n, n * n);
      m.setFromTriplets(sparse_list.begin(), sparse_list.end());
      return m;
    }

    class cholesky_decompose_v_vari : public vari {
    public:
      int M_;  // A.rows() = A.cols() = B.rows()
      double* L_;
      vari** _variRefA;
      vari** _variRefL;
      vari** _dummy;

      cholesky_decompose_v_vari(const Eigen::Matrix<var, -1, -1>& A,
                                const Eigen::Matrix<double, -1, -1>& L_A)
        : vari(0.0),
          M_(A.rows()),
          L_(reinterpret_cast<double*>
             (stan::math::ChainableStack::memalloc_
              .alloc(sizeof(double) * A.rows() * (A.rows() + 1) / 2))),
          _variRefA(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*) * A.rows() * (A.cols() + 1) / 2))),
          _variRefL(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*) * A.rows() * (A.rows() + 1) / 2))),
          _dummy(reinterpret_cast<vari**>
                    (stan::math::ChainableStack::memalloc_
                     .alloc(sizeof(vari*))))  {

        size_t pos = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            _variRefA[pos++] = A(i, j).vi_;
          }
        }

        pos = 0;
        for (size_type j = 0; j < M_; ++j) {
          for (size_type i = j; i < M_; ++i) {
            L_[pos] = L_A(i, j);
            _variRefL[pos] = new vari(L_[pos], false);
            pos++;
          }
        }
        _dummy[0] = new vari(0.0, false);
      }

      virtual void chain() {
        using Eigen::Matrix;
        using Eigen::Dynamic;
        using Eigen::Lower;
        using Eigen::OnTheRight;
        using Eigen::Map;
        Eigen::Matrix<double, 1, Dynamic> adjL(M_ * (M_ + 1) / 2);
        Eigen::Matrix<double, Dynamic, Dynamic> LA(M_, M_);

        size_t pos = -1;
        for (size_type j = 0; j < M_; ++j)
          for (size_type i = j; i < M_; ++i) {
            ++pos;
            adjL(pos) = _variRefL[pos]->adj_;
            LA(i, j) = L_[pos];
          }

        Eigen::SparseMatrix<double> dAdL = id_p_commutation_matrix(M_) * chol_kron_eye(LA);
        Eigen::SparseMatrix<double> elim = elimination_matrix<double>(M_);
        Eigen::SparseMatrix<double> res = elim * dAdL * Eigen::SparseMatrix<double>(elim.transpose());
        Eigen::MatrixXd dvAdvL = Eigen::MatrixXd(res);
        dvAdvL.triangularView<Lower>().solveInPlace<Eigen::OnTheRight>(adjL);

        for (size_type j = 0; j < M_ * (M_ + 1) / 2; ++j)
          _variRefA[j]->adj_ += adjL(j);
      }
    };

    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    cholesky_decompose(const Eigen::Matrix<var, -1, -1> &A) {
      Eigen::Matrix<var, -1, -1> L(A.rows(), A.cols());
      Eigen::Matrix<double, -1, -1> L_A(A.rows(), A.cols());
      
      L_A = cholesky_decompose(value_of_rec(A));
      // NOTE: this is not a memory leak, this vari is used in the
      // expression graph to evaluate the adjoint, but is not needed
      // for the returned matrix.  Memory will be cleaned up with the
      // arena allocator.
      cholesky_decompose_v_vari *baseVari
        = new cholesky_decompose_v_vari(A, L_A);

      size_t pos = 0;
      for (size_type j = 0; j < L.cols(); ++j) {
        for (size_type i = j; i < L.rows(); ++i) 
          L(i, j).vi_ = baseVari->_variRefL[pos++];
        for (size_type k = 0; k < j; ++k) 
          L(k, j).vi_ = baseVari->_dummy[0];
      }
      return L;
    }

  }
}
#endif
