#ifndef STAN_MATH_REV_MAT_FUN_QUAD_FORM_HPP
#define STAN_MATH_REV_MAT_FUN_QUAD_FORM_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/quad_form.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>

namespace stan {
namespace math {

namespace {
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class quad_form_vari : public vari {
 protected:
  inline void chainA(double* A_mem, const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {}
  inline void chainB(double* B_mem, const Eigen::Matrix<double, Ra, Ca>& Ad,
                     const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {}

  inline void chainA(var* A_mem, const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {
    Eigen::Matrix<double, Ra, Ca> adjA(Bd * adjC * Bd.transpose());
    for (int j = 0; j < A_rows_; j++) {
      for (int i = 0; i < A_rows_; i++) {
        A_mem[i + j * A_rows_].vi_->adj_ += adjA(i, j);
      }
    }
  }
  inline void chainB(var* B_mem, const Eigen::Matrix<double, Ra, Ca>& Ad,
                     const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {
    Eigen::Matrix<double, Ra, Ca> adjB(Ad * Bd * adjC.transpose()
                                       + Ad.transpose() * Bd * adjC);
    for (int j = 0; j < B_cols_; j++)
      for (int i = 0; i < B_rows_; i++)
        B_mem[i + j * B_rows_].vi_->adj_ += adjB(i, j);
  }

  inline void chainAB(const Eigen::Matrix<double, Ra, Ca>& Ad,
                      const Eigen::Matrix<double, Rb, Cb>& Bd,
                      const Eigen::Matrix<double, Cb, Cb>& adjC) {
    chainA(A_mem_, Bd, adjC);
    chainB(B_mem_, Ad, Bd, adjC);
  }

 public:
  quad_form_vari(const Eigen::Matrix<Ta, Ra, Ca>& A,
                 const Eigen::Matrix<Tb, Rb, Cb>& B, bool symmetric = false)
      : vari(0.0),
        A_rows_(A.rows()),
        B_rows_(B.rows()),
        B_cols_(B.cols()),
        C_rows_(B.cols()),
        A_mem_(ChainableStack::instance().memalloc_.alloc_array<Ta>(A.size())),
        B_mem_(ChainableStack::instance().memalloc_.alloc_array<Tb>(B.size())),
        C_mem_(ChainableStack::instance().memalloc_.alloc_array<var>(
            C_rows_ * C_rows_)) {
    for (int i = 0; i < A.size(); ++i)
      A_mem_[i] = A.data()[i];
    for (int i = 0; i < B.size(); ++i)
      B_mem_[i] = B.data()[i];

    const Eigen::Matrix<double, Ra, Ca>& Ad = value_of(A);
    const Eigen::Matrix<double, Ra, Ca>& Bd = value_of(B);
    Eigen::Matrix<double, Cb, Cb> Cd(Bd.transpose() * Ad * Bd);

    for (int j = 0; j < C_rows_; ++j) {
      for (int i = 0; i < C_rows_; ++i) {
        if (symmetric) {
          C_mem_[i + C_rows_ * j]
              = var(new vari(0.5 * (Cd(i, j) + Cd(j, i)), false));
        } else {
          C_mem_[i + C_rows_ * j] = var(new vari(Cd(i, j), false));
        }
      }
    }
  }

  virtual void chain() {
    Eigen::Matrix<double, Ra, Ca> Ad(A_rows_, A_rows_);
    Eigen::Matrix<double, Rb, Cb> Bd(B_rows_, B_cols_);
    Eigen::Map<Eigen::Matrix<var, Cb, Cb> > C(C_mem_, C_rows_, C_rows_);

    for (int j = 0; j < A_rows_; ++j)
      for (int i = 0; i < A_rows_; ++i)
        Ad(i, j) = value_of(A_mem_[i + A_rows_ * j]);

    for (int j = 0; j < B_cols_; ++j)
      for (int i = 0; i < B_rows_; ++i)
        Bd(i, j) = value_of(B_mem_[i + B_rows_ * j]);

    Eigen::Matrix<double, Cb, Cb> adjC(C_rows_, C_rows_);
    for (int j = 0; j < C_rows_; ++j)
      for (int i = 0; i < C_rows_; ++i)
        adjC(i, j) = C(i, j).vi_->adj_;

    chainAB(Ad, Bd, adjC);
  }

  int A_rows_;  // A is square
  int B_rows_;
  int B_cols_;
  int C_rows_;  // C is square
  Ta* A_mem_;
  Tb* B_mem_;
  var* C_mem_;
};
}  // namespace

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
inline typename boost::enable_if_c<boost::is_same<Ta, var>::value
                                       || boost::is_same<Tb, var>::value,
                                   Eigen::Matrix<var, Cb, Cb> >::type
quad_form(const Eigen::Matrix<Ta, Ra, Ca>& A,
          const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>* baseVari
      = new quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>(A, B);

  Eigen::Matrix<var, Cb, Cb> C(baseVari->C_rows_, baseVari->C_rows_);
  for (int i = 0; i < baseVari->C_rows_ * baseVari->C_rows_; ++i)
    C.data()[i] = baseVari->C_mem_[i];

  return C;
}

template <typename Ta, int Ra, int Ca, typename Tb, int Rb>
inline typename boost::enable_if_c<
    boost::is_same<Ta, var>::value || boost::is_same<Tb, var>::value, var>::type
quad_form(const Eigen::Matrix<Ta, Ra, Ca>& A,
          const Eigen::Matrix<Tb, Rb, 1>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>* baseVari
      = new quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>(A, B);

  return baseVari->C_mem_[0];
}

}  // namespace math
}  // namespace stan
#endif
