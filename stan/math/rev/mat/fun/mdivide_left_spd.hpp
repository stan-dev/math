#ifndef STAN_MATH_REV_MAT_FUN_MDIVIDE_LEFT_SPD_HPP
#define STAN_MATH_REV_MAT_FUN_MDIVIDE_LEFT_SPD_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <vector>

namespace stan {
namespace math {

namespace {

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class mdivide_left_spd_vari : public vari {
 protected:
  Ta* A_mem_;
  const double* Ad_llt_mem_;
  int A_rows_;  // A is square
  Tb* b_mem_;
  int b_rows_;
  int b_cols_;
  int C_rows_;
  int C_cols_;
  vari** C_mem_;
  const double* Cd_mem_;

  static inline void chainA(
      const double* A_mem, Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Cd,
      const Eigen::Matrix<double, Rb, Cb>& adjb) {}

  static inline void chainb(const double* b_mem,
                            const Eigen::Matrix<double, Rb, Cb>& adjb) {}

  static inline void chainA(
      vari** A_mem, Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Cd,
      const Eigen::Matrix<double, Rb, Cb>& adjb) {
    Eigen::Matrix<double, Ra, Ca> adjA = -adjb * Cd.transpose();

    for (int i = 0; i < adjA.size(); ++i) {
      A_mem[i]->adj_ += adjA(i);
    }
  }

  static inline void chainb(vari** b_mem,
                            const Eigen::Matrix<double, Rb, Cb>& adjb) {
    for (int i = 0; i < adjb.size(); ++i) {
      b_mem[i]->adj_ += adjb(i);
    }
  }

 public:
  mdivide_left_spd_vari(const double* Cd_mem, Ta* A_mem,
                        const double* Ad_llt_mem, int A_rows, Tb* b_mem,
                        int b_rows, int b_cols)
      : vari(Cd_mem[0]),
        A_mem_(A_mem),
        Ad_llt_mem_(Ad_llt_mem),
        A_rows_(A_rows),
        b_mem_(b_mem),
        b_rows_(b_rows),
        b_cols_(b_cols),
        C_rows_(b_rows),
        C_cols_(b_cols),
        C_mem_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            A_rows * b_cols)),
        Cd_mem_(Cd_mem) {
    C_mem_[0] = this;
    for (int i = 1; i < C_rows_ * C_cols_; ++i) {
      C_mem_[i] = new vari(Cd_mem[i], false);
    }
  }

  Eigen::Matrix<var, Rb, Cb> get_output() {
    Eigen::Matrix<var, Rb, Cb> output(b_rows_, b_cols_);
    for (int i = 0; i < b_rows_ * b_cols_; ++i)
      output(i) = C_mem_[i];
    return output;
  }

  virtual void chain() {
    Eigen::Map<const Eigen::Matrix<double, Ra, Ca> > Ad_llt(Ad_llt_mem_,
                                                            A_rows_, A_rows_);
    Eigen::Map<const Eigen::Matrix<double, Rb, Cb> > Cd(Cd_mem_, C_rows_,
                                                        C_cols_);

    Eigen::Matrix<double, Rb, Cb> adjb(b_rows_, b_cols_);
    for (int i = 0; i < C_cols_ * C_rows_; ++i) {
      adjb(i) = C_mem_[i]->adj_;
    }

    Ad_llt.template triangularView<Eigen::Lower>().solveInPlace(adjb);
    Ad_llt.template triangularView<Eigen::Upper>().solveInPlace(adjb);

    chainA(A_mem_, Cd, adjb);
    chainb(b_mem_, adjb);
  }
};

template <int Ra, int Ca>
const Eigen::Matrix<double, Ra, Ca>& get_Ad(
    const Eigen::Matrix<double, Ra, Ca>& A) {
  return A;
}

template <int Ra, int Ca>
Eigen::Matrix<double, Ra, Ca> get_Ad(const Eigen::Matrix<var, Ra, Ca>& A) {
  Eigen::Matrix<double, Ra, Ca> Ad(A.rows(), A.cols());

  for (int i = 0; i < A.size(); ++i)
    Ad(i) = A(i).val();

  return Ad;
}

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
Eigen::Matrix<var, Ca, Cb> mdivide_left_spd_rev(
    const Eigen::Matrix<Ta, Ra, Ca>& A, const Eigen::Matrix<Tb, Rb, Cb>& b) {
  check_square("mdivide_left_spd", "A", A);
  check_multiplicable("mdivide_left_spd", "A", A, "b", b);

  auto A_mem = build_vari_pointer_array_if_necessary(A.data(), A.size());
  auto b_mem = build_vari_pointer_array_if_necessary(b.data(), b.size());

  auto Ad_llt = get_Ad(A).llt();
  double* Cd_mem = build_double_array(b_mem, b.size());

  Eigen::Map<Eigen::Matrix<double, Rb, Cb> > Cd(Cd_mem, b.rows(), b.cols());

  Ad_llt.solveInPlace(Cd);

  double* Ad_llt_mem
      = ChainableStack::instance().memalloc_.alloc_array<double>(A.size());
  Eigen::Map<Eigen::Matrix<double, Ra, Ca> > Ad_llt_matrix(Ad_llt_mem, A.rows(),
                                                           A.cols());

  Ad_llt_matrix.template triangularView<Eigen::Lower>() = Ad_llt.matrixL();
  Ad_llt_matrix.template triangularView<Eigen::Upper>()
      = Ad_llt_matrix.transpose();

  auto baseVari
      = new mdivide_left_spd_vari<std::remove_pointer_t<decltype(A_mem)>, Ra,
                                  Ca, std::remove_pointer_t<decltype(b_mem)>,
                                  Rb, Cb>(Cd_mem, A_mem, Ad_llt_mem, A.rows(),
                                          b_mem, b.rows(), b.cols());

  return baseVari->get_output();
}
}  // namespace

template <int Ra, int Ca, int Rb, int Cb>
Eigen::Matrix<var, Ca, Cb> mdivide_left_spd(
    const Eigen::Matrix<var, Ra, Ca>& A, const Eigen::Matrix<var, Rb, Cb>& b) {
  return mdivide_left_spd_rev(A, b);
}

template <int Ra, int Ca, int Rb, int Cb>
Eigen::Matrix<var, Ca, Cb> mdivide_left_spd(
    const Eigen::Matrix<var, Ra, Ca>& A,
    const Eigen::Matrix<double, Rb, Cb>& b) {
  return mdivide_left_spd_rev(A, b);
}

template <int Ra, int Ca, int Rb, int Cb>
Eigen::Matrix<var, Ca, Cb> mdivide_left_spd(
    const Eigen::Matrix<double, Ra, Ca>& A,
    const Eigen::Matrix<var, Rb, Cb>& b) {
  return mdivide_left_spd_rev(A, b);
}
}  // namespace math
}  // namespace stan
#endif
