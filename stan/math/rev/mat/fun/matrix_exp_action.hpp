#ifndef STAN_MATH_REV_MAT_FUN_MATRIX_EXP_ACTION_HPP
#define STAN_MATH_REV_MAT_FUN_MATRIX_EXP_ACTION_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/matrix_exp_action_handler.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>

namespace stan {
namespace math {

  Eigen::MatrixXd exp_action_chain_dv(double* Ad,
                                      const int n,
                                      const Eigen::MatrixXd& adjexpAB,
                                      const double t) {
      using Eigen::Map;
      matrix_exp_action_handler handle;
      return handle.action(Map<Eigen::MatrixXd>(Ad, n, n).transpose(), 
                             adjexpAB, t);
  }

  Eigen::MatrixXd exp_action_chain_vd(double* Ad,
                                      const int n,
                                      const Eigen::MatrixXd& adjexpAB,
                                      const double t) {
    Eigen::MatrixXd adjA = Eigen::MatrixXd::Zero(n, n);
    // TODO(yizhang): a better way like complex step approximation
    try {
      start_nested();
      Eigen::Matrix<stan::math::var,
                    Eigen::Dynamic,
                    Eigen::Dynamic> A_temp(n, n);
      std::vector<stan::math::var> Avar;
      std::vector<double> g;
      for (int i = 0; i < n * n; ++i) {
        A_temp(i) = to_var(Ad[i]);
        Avar.push_back(A_temp(i));
      }
      Eigen::Matrix<stan::math::var,
                    Eigen::Dynamic,
                    Eigen::Dynamic> expA = matrix_exp(A_temp);
      for (size_type i = 0; i < expA.size(); ++i) {
        stan::math::set_zero_all_adjoints_nested();
        expA.coeffRef(i).grad(Avar, g);
        for (size_type j = 0; j < n * n; ++j) {
          adjA(j) += adjexpAB(i) * g[j];
        }
      }
    } catch (const std::exception& e) {
      recover_memory_nested();
      throw;
    }
    recover_memory_nested();
    return adjA;
}

  template <typename Ta, int N, typename Tb, int Cb>
  class matrix_exp_action_vari : public vari {
  public:
    int n_;
    int B_cols_;
    int A_size_;
    int B_size_;
    double t_;
    double* Ad_;
    double* Bd_;
    vari** variRefA_;
    vari** variRefB_;
    vari** variRefexpAB_;
    
    matrix_exp_action_vari(const Eigen::Matrix<Ta, N, N>& A,
                           const Eigen::Matrix<Tb, N, Cb>& B,
                           const double& t)
      : vari(0.0),
        n_(A.rows()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        t_(t),
        Ad_(ChainableStack::memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::memalloc_.alloc_array<double>(B_size_)),
        variRefA_(ChainableStack::memalloc_.alloc_array<vari*>(A_size_)),
        variRefB_(ChainableStack::memalloc_.alloc_array<vari*>(B_size_)),
        variRefexpAB_(
            ChainableStack::memalloc_.alloc_array<vari*>(B_size_)) {
      using Eigen::Map;
      using Eigen::MatrixXd;
      for (size_type i = 0; i < A.size(); ++i) {
        variRefA_[i] = A.coeffRef(i).vi_;
        Ad_[i] = A.coeffRef(i).val();
      }
      for (size_type i = 0; i < B.size(); ++i) {
        variRefB_[i] = B.coeffRef(i).vi_;
        Bd_[i] = B.coeffRef(i).val();
      }
      matrix_exp_action_handler handle;
      MatrixXd expAB = handle.action(Map<MatrixXd>(Ad_, n_, n_),
                                     Map<MatrixXd>(Bd_, n_, B_cols_),
                                     t_);
      for (size_type i = 0; i < expAB.size(); ++i)
        variRefexpAB_[i] = new vari(expAB.coeffRef(i), false);
    }

    virtual void chain() {
      using Eigen::Map;
      using Eigen::MatrixXd;
      MatrixXd adjexpAB(n_, B_cols_);

      for (size_type i = 0; i < adjexpAB.size(); ++i)
        adjexpAB(i) = variRefexpAB_[i]->adj_;

      MatrixXd adjA = exp_action_chain_dv(Ad_, n_, adjexpAB, t_);
      MatrixXd adjB = exp_action_chain_dv(Ad_, n_, adjexpAB, t_);

      for (size_type i = 0; i < A_size_; ++i) {
        variRefA_[i]->adj_ += adjA(i);
      }
      for (size_type i = 0; i < B_size_; ++i) {
        variRefB_[i]->adj_ += adjB(i);
      }
    }
  };
  

  template <int N, typename Tb, int Cb>
  class matrix_exp_action_vari<double, N, Tb, Cb> : public vari {
  public:
    int n_;
    int B_cols_;
    int A_size_;
    int B_size_;
    double t_;
    double* Ad_;
    double* Bd_;
    vari** variRefB_;
    vari** variRefexpAB_;

    matrix_exp_action_vari(const Eigen::Matrix<double, N, N>& A,
                           const Eigen::Matrix<Tb, N, Cb>& B,
                           const double& t)
      : vari(0.0),
        n_(A.rows()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        t_(t),
        Ad_(ChainableStack::memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::memalloc_.alloc_array<double>(B_size_)),
        variRefB_(ChainableStack::memalloc_.alloc_array<vari*>(B_size_)),
        variRefexpAB_(
                      ChainableStack::memalloc_.alloc_array<vari*>(B_size_)) {
      using Eigen::Map;
      using Eigen::MatrixXd;
      for (size_type i = 0; i < A.size(); ++i)
        Ad_[i] = A.coeffRef(i);
      for (size_type i = 0; i < B.size(); ++i) {
        variRefB_[i] = B.coeffRef(i).vi_;
        Bd_[i] = B.coeffRef(i).val();
      }
      matrix_exp_action_handler handle;
      MatrixXd expAB = handle.action(Map<MatrixXd>(Ad_, n_, n_),
                                     Map<MatrixXd>(Bd_, n_, B_cols_),
                                     t_);
      for (size_type i = 0; i < expAB.size(); ++i)
        variRefexpAB_[i] = new vari(expAB.coeffRef(i), false);
    }

    virtual void chain() {
      using Eigen::Map;
      using Eigen::MatrixXd;
      MatrixXd adjexpAB(n_, B_cols_);
      // MatrixXd adjB(n_, B_cols_);

      for (size_type i = 0; i < adjexpAB.size(); ++i)
        adjexpAB(i) = variRefexpAB_[i]->adj_;

      MatrixXd adjB = exp_action_chain_dv(Ad_, n_, adjexpAB, t_);

      for (size_type i = 0; i < B_size_; ++i) {
        variRefB_[i]->adj_ += adjB(i);
      }
    }
  };
  
  template <typename Ta, int N, int Cb>
  class matrix_exp_action_vari<Ta, N, double, Cb> : public vari {
  public:
    int n_;
    int B_cols_;
    int A_size_;
    int B_size_;
    double t_;
    double* Ad_;
    double* Bd_;
    vari** variRefA_;
    vari** variRefexpAB_;

    matrix_exp_action_vari(const Eigen::Matrix<Ta, N, N>& A,
                           const Eigen::Matrix<double, N, Cb>& B,
                           const double& t)
      : vari(0.0),
        n_(A.rows()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        t_(t),
        Ad_(ChainableStack::memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::memalloc_.alloc_array<double>(B_size_)),
        variRefA_(ChainableStack::memalloc_.alloc_array<vari*>(A_size_)),
        variRefexpAB_(
                      ChainableStack::memalloc_.alloc_array<vari*>(B_size_)) {
      using Eigen::Map;
      using Eigen::MatrixXd;
      for (size_type i = 0; i < A.size(); ++i) {
        variRefA_[i] = A.coeffRef(i).vi_;
        Ad_[i] = A.coeffRef(i).val();
      }
      for (size_type i = 0; i < B.size(); ++i) {
        Bd_[i] = B.coeffRef(i);
      }
      matrix_exp_action_handler handle;
      MatrixXd expAB = handle.action(Map<MatrixXd>(Ad_, n_, n_),
                                     Map<MatrixXd>(Bd_, n_, B_cols_),
                                     t_);
      for (size_type i = 0; i < expAB.size(); ++i)
        variRefexpAB_[i] = new vari(expAB.coeffRef(i), false);

    }

    virtual void chain() {
      using Eigen::Map;
      using Eigen::MatrixXd;
      MatrixXd adjexpAB(n_, B_cols_);

      for (size_type i = 0; i < adjexpAB.size(); ++i)
        adjexpAB(i) = variRefexpAB_[i]->adj_;

      MatrixXd adjA = exp_action_chain_dv(Ad_, n_, adjexpAB, t_);

      for (size_type i = 0; i < A_size_; ++i) {
        variRefA_[i]->adj_ += adjA(i);
      }
    }
  };

  template <typename Ta, int N, typename Tb, int Cb>
  inline
  typename boost::enable_if_c<boost::is_same<Ta, var>::value
                              || boost::is_same<Tb, var>::value,
                              Eigen::Matrix<var, N, Cb> >::type
  matrix_exp_action(const Eigen::Matrix<Ta, N, N>& A,
                    const Eigen::Matrix<Tb, N, Cb>& B,
                    const double& t = 1.0) {
    const char* fun = "matrix_exp_action";
    check_multiplicable(fun, "A", A, "B", B);
    // check_not_nan(fun, "A", A);
    // check_not_nan(fun, "B", B);

    matrix_exp_action_vari<Ta, N, Tb, Cb>* baseVari
      = new matrix_exp_action_vari<Ta, N, Tb, Cb>(A, B, t);
    Eigen::Matrix<var, N, Cb> expAB_v(A.rows(), B.cols());
    for (size_type i = 0; i < expAB_v.size(); ++i) {
      expAB_v.coeffRef(i).vi_ = baseVari->variRefexpAB_[i];
    }
    return expAB_v;
  }

  template <int N, int Cb>
  inline
  Eigen::Matrix<double, N, Cb> 
  matrix_exp_action(const Eigen::Matrix<double, N, N>& A,
                    const Eigen::Matrix<double, N, Cb>& B,
                    const double& t = 1.0) {
    const char* fun = "matrix_exp_action";
    check_multiplicable(fun, "A", A, "B", B);
    // check_not_nan(fun, "A", A);
    // check_not_nan(fun, "B", B);

    Eigen::Matrix<double, N, Cb> expAB;
    matrix_exp_action_handler handle;
    expAB = handle.action(A, B, t);
    return expAB;
  }

}  // namespace math
}  // namespace stan

#endif
