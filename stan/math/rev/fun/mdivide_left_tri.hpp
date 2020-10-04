#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_TRI_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_TRI_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl.hpp>
#endif

namespace stan {
namespace math {

namespace internal {

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_vv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefA_;
  vari **variRefB_;
  vari **variRefC_;

  mdivide_left_tri_vv_vari(const Eigen::Matrix<var, R1, C1> &A,
                           const Eigen::Matrix<var, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        C_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * B.rows()
                                                       * B.cols()))),
        variRefA_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * (A.rows() + 1) / 2))),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;

    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = j; i < M_; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    }

    Map<matrix_d> c_map(C_, M_, N_);
    Map<matrix_d> a_map(A_, M_, M_);
    a_map = A.val();
    c_map = B.val();
    Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
#ifdef STAN_OPENCL
    if (A.rows()
        >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
      matrix_cl<double> A_cl(a_map, from_eigen_uplo_type(TriView));
      matrix_cl<double> C_cl(c_map);
      C_cl = tri_inverse(A_cl) * C_cl;
      c_map = from_matrix_cl(C_cl);
    } else {
#endif
      c_map = a_map.template triangularView<TriView>().solve(c_map);
#ifdef STAN_OPENCL
    }
#endif

    Map<matrix_vi>(variRefC_, M_, N_)
        = c_map.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    matrix_d adjA;
    matrix_d adjB;

#ifdef STAN_OPENCL
    if (M_ >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
      matrix_cl<double> A_cl(A_, M_, M_, from_eigen_uplo_type(TriView));
      matrix_cl<double> C_cl(C_, M_, N_);
      matrix_cl<double> variRefC_cl(Map<matrix_vi>(variRefC_, M_, N_).adj());
      matrix_cl<double> adjB_cl = transpose(tri_inverse(A_cl)) * variRefC_cl;
      matrix_cl<double> adjA_cl = adjB_cl * transpose(C_cl) * -1.0;
      adjA = from_matrix_cl(adjA_cl);
      adjB = from_matrix_cl(adjB_cl);
    } else {
#endif
      adjB = Map<matrix_d>(A_, M_, M_)
                 .template triangularView<TriView>()
                 .transpose()
                 .solve(Map<matrix_vi>(variRefC_, M_, N_).adj());
      adjA = -adjB * Map<matrix_d>(C_, M_, N_).transpose();
#ifdef STAN_OPENCL
    }
#endif
    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = j; i < adjA.rows(); i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    }
    Map<matrix_vi>(variRefB_, M_, N_).adj() += adjB;
  }
};

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_dv_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefB_;
  vari **variRefC_;

  mdivide_left_tri_dv_vari(const Eigen::Matrix<double, R1, C1> &A,
                           const Eigen::Matrix<var, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        C_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * B.rows()
                                                       * B.cols()))),
        variRefB_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;

    Map<matrix_d>(A_, M_, M_) = A;
    Map<matrix_vi>(variRefB_, M_, N_) = B.vi();
    Map<matrix_d> c_map(C_, M_, N_);
    c_map = B.val();

#ifdef STAN_OPENCL
    if (A.rows()
        >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
      matrix_cl<double> A_cl(A, from_eigen_uplo_type(TriView));
      matrix_cl<double> C_cl(c_map);
      C_cl = tri_inverse(A_cl) * C_cl;
      c_map = from_matrix_cl(C_cl);
    } else {
#endif
      c_map = Map<matrix_d>(A_, M_, M_)
                  .template triangularView<TriView>()
                  .solve(c_map);
#ifdef STAN_OPENCL
    }
#endif
    Map<matrix_vi>(variRefC_, M_, N_)
        = c_map.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
#ifdef STAN_OPENCL
    if (M_ >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
      matrix_cl<double> A_cl(A_, M_, M_, from_eigen_uplo_type(TriView));
      matrix_cl<double> C_cl(Map<matrix_vi>(variRefC_, M_, N_).adj());
      A_cl = transpose(tri_inverse(A_cl));
      matrix_cl<double> res_cl = A_cl * C_cl;
      Map<matrix_vi>(variRefB_, M_, N_).adj() += from_matrix_cl(res_cl);
    } else {
#endif
      Map<matrix_vi>(variRefB_, M_, N_).adj()
          += Map<matrix_d>(A_, M_, M_)
                 .template triangularView<TriView>()
                 .transpose()
                 .solve(Map<matrix_vi>(variRefC_, M_, N_).adj());
#ifdef STAN_OPENCL
    }
#endif
  }
};

template <Eigen::UpLoType TriView, int R1, int C1, int R2, int C2>
class mdivide_left_tri_vd_vari : public vari {
 public:
  int M_;  // A.rows() = A.cols() = B.rows()
  int N_;  // B.cols()
  double *A_;
  double *C_;
  vari **variRefA_;
  vari **variRefC_;

  mdivide_left_tri_vd_vari(const Eigen::Matrix<var, R1, C1> &A,
                           const Eigen::Matrix<double, R2, C2> &B)
      : vari(0.0),
        M_(A.rows()),
        N_(B.cols()),
        A_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * A.rows()
                                                       * A.cols()))),
        C_(reinterpret_cast<double *>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(double) * B.rows()
                                                       * B.cols()))),
        variRefA_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * A.rows()
                                                       * (A.rows() + 1) / 2))),
        variRefC_(reinterpret_cast<vari **>(
            ChainableStack::instance_->memalloc_.alloc(sizeof(vari *) * B.rows()
                                                       * B.cols()))) {
    using Eigen::Map;
    using Eigen::Matrix;

    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = j; i < M_; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < M_; j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++] = A(i, j).vi_;
        }
      }
    }
    Map<matrix_d> Ad(A_, M_, M_);
    Map<matrix_d> Cd(C_, M_, N_);
    Ad = A.val();
#ifdef STAN_OPENCL
    if (M_ >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
      matrix_cl<double> A_cl(Ad, from_eigen_uplo_type(TriView));
      matrix_cl<double> B_cl(B);
      B_cl = tri_inverse(A_cl) * B_cl;
      Cd = from_matrix_cl(B_cl);
    } else {
#endif
      Cd = Ad.template triangularView<TriView>().solve(B);
#ifdef STAN_OPENCL
    }
#endif
    Map<matrix_vi>(variRefC_, M_, N_)
        = Cd.unaryExpr([](double x) { return new vari(x, false); });
  }

  virtual void chain() {
    using Eigen::Map;
    using Eigen::Matrix;
    Matrix<double, R1, C1> adjA(M_, M_);
    Matrix<double, R1, C2> adjC(M_, N_);

    adjC = Map<matrix_vi>(variRefC_, M_, N_).adj();
#ifdef STAN_OPENCL
    if (M_ >= opencl_context.tuning_opts().tri_inverse_size_worth_transfer) {
      matrix_cl<double> A_cl(A_, M_, M_, from_eigen_uplo_type(TriView));
      matrix_cl<double> C_cl(C_, M_, N_);
      matrix_cl<double> adjC_cl(adjC);
      A_cl = transpose(tri_inverse(A_cl));
      matrix_cl<double> adjA_cl = A_cl * (adjC_cl * transpose(C_cl)) * -1.0;
      adjA = from_matrix_cl(adjA_cl);
    } else {
#endif
      adjA.noalias()
          = -Map<Matrix<double, R1, C1>>(A_, M_, M_)
                 .template triangularView<TriView>()
                 .transpose()
                 .solve(adjC
                        * Map<Matrix<double, R1, C2>>(C_, M_, N_).transpose());
#ifdef STAN_OPENCL
    }
#endif
    size_t pos = 0;
    if (TriView == Eigen::Lower) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = j; i < adjA.rows(); i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    } else if (TriView == Eigen::Upper) {
      for (size_type j = 0; j < adjA.cols(); j++) {
        for (size_type i = 0; i < j + 1; i++) {
          variRefA_[pos++]->adj_ += adjA(i, j);
        }
      }
    }
  }
};
}  // namespace internal

template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_eigen_t<T1, T2> * = nullptr,
          require_any_vt_var<T1, T2> * = nullptr>
inline auto mdivide_left_tri(const T1 &A, const T2 &B) {
  using ret_type = promote_scalar_t<var, decltype(A * B)>;

  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "B", B);

  if (A.size() == 0) {
    return ret_type(0, B.cols());
  }

  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);

  arena_matrix<promote_scalar_t<var, T1>> arena_A;
  arena_matrix<promote_scalar_t<var, T2>> arena_B;

  if (!is_constant<T1>::value) {
    arena_A = A_ref;
  }

  if (!is_constant<T2>::value) {
    arena_B = B_ref;
  }

  arena_t<promote_scalar_t<double, T1>> arena_A_val
      = value_of(A_ref).template triangularView<TriView>();
  arena_t<ret_type> res
      = arena_A_val.template triangularView<TriView>().solve(value_of(B_ref));

  reverse_pass_callback([arena_A, arena_B, arena_A_val, res]() mutable {
    promote_scalar_t<double, T2> adjB
        = arena_A_val.template triangularView<TriView>().transpose().solve(
            res.adj());

    if (!is_constant<T2>::value)
      arena_B.adj() += adjB;

    if (!is_constant<T1>::value)
      arena_A.adj() -= (adjB * res.val().transpose().eval())
	.template triangularView<TriView>();
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
