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
  Ta* A_mem_;
  const double* Ad_mem_;
  int A_rows_;  // A is square
  Tb* B_mem_;
  const double* Bd_mem_;
  int B_rows_;
  int B_cols_;
  int C_rows_;  // C is square
  vari** C_mem_;

  inline void chainA(const double* A_mem,
                     Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {}
  inline void chainB(const double* B_mem,
                     Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >& Ad,
                     Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {}

  inline void chainA(vari** A_mem,
                     Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {
    Eigen::Matrix<double, Ra, Ca> adjA(Bd * adjC * Bd.transpose());
    for (int i = 0; i < adjA.size(); ++i) {
      A_mem[i]->adj_ += adjA.data()[i];
    }
  }
  inline void chainB(vari** B_mem,
                     Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >& Ad,
                     Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {
    Eigen::Matrix<double, Ra, Ca> adjB(Ad * Bd * adjC.transpose()
                                       + Ad.transpose() * Bd * adjC);
    for (int i = 0; i < adjB.size(); ++i) {
      B_mem[i]->adj_ += adjB.data()[i];
    }
  }

 public:
  quad_form_vari(const Eigen::Matrix<double, Ca, Cb>& Cd, Ta* A_mem,
                 const double* Ad_mem, int A_rows, Tb* B_mem,
                 const double* Bd_mem, int B_rows, int B_cols,
                 bool symmetric = false)
      : vari(Cd(0, 0)),
        A_mem_(A_mem),
        Ad_mem_(Ad_mem),
        A_rows_(A_rows),
        B_mem_(B_mem),
        Bd_mem_(Bd_mem),
        B_rows_(B_rows),
        B_cols_(B_cols),
        C_rows_(B_cols),
        C_mem_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
            C_rows_ * C_rows_)) {
    for (int j = 0; j < C_rows_; ++j) {
      for (int i = 0; i < C_rows_; ++i) {
        if (j == 0 && i == 0) {
          C_mem_[0] = this;
        } else {
          if (symmetric) {
            C_mem_[i + C_rows_ * j]
                = new vari(0.5 * (Cd(i, j) + Cd(j, i)), false);
          } else {
            C_mem_[i + C_rows_ * j] = new vari(Cd(i, j), false);
          }
        }
      }
    }
  }

  Eigen::Matrix<var, Cb, Cb> getOutput() {
    Eigen::Matrix<var, Cb, Cb> C(C_rows_, C_rows_);
    for (int i = 0; i < C_rows_ * C_rows_; ++i)
      C.data()[i] = var(C_mem_[i]);
    return C;
  }

  virtual void chain() {
    Eigen::Map<const Eigen::Matrix<double, Ra, Ca> > Ad(Ad_mem_, A_rows_,
                                                        A_rows_);
    Eigen::Map<const Eigen::Matrix<double, Rb, Cb> > Bd(Bd_mem_, B_rows_,
                                                        B_cols_);

    Eigen::Matrix<double, Cb, Cb> adjC(C_rows_, C_rows_);
    for (int i = 0; i < C_rows_ * C_rows_; ++i)
      adjC.data()[i] = C_mem_[i]->adj_;

    chainA(A_mem_, Bd, adjC);
    chainB(B_mem_, Ad, Bd, adjC);
  }
};

/**
 * quad_form_inner computes the quadratic form B^T * A * B
 *
 * The scalars in A and B can be either doubles or vars. The scalars in
 * at least one of them should be vars.
 *
 * @tparam Ta Type of scalar in A
 * @tparam Ra Eigen row type of A
 * @tparam Ca Eigen column type of A
 * @tparam Tb Type of scalar in B
 * @tparam Rb Eigen row type of B
 * @tparam Cb Eigen column type of B
 * @param A A matrix in quadratic form
 * @param B B matrix (or vector) in quadratic form
 * @return an Eigen Matrix of vars containing the quadratic form (B^T * A * B)
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
Eigen::Matrix<var, Cb, Cb> quad_form_inner(const Eigen::Matrix<Ta, Ra, Ca>& A,
                                           const Eigen::Matrix<Tb, Rb, Cb>& B,
                                           bool symmetric) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  auto A_mem = build_vari_pointer_array_if_necessary(A.data(), A.size());
  auto B_mem = build_vari_pointer_array_if_necessary(B.data(), B.size());

  const double* Ad_mem = build_double_array(A_mem, A.size());
  const double* Bd_mem = build_double_array(B_mem, B.size());

  Eigen::Map<const Eigen::Matrix<double, Ra, Ca> > Ad(Ad_mem, A.rows(),
                                                      A.cols());
  Eigen::Map<const Eigen::Matrix<double, Rb, Cb> > Bd(Bd_mem, B.rows(),
                                                      B.cols());
  Eigen::Matrix<double, Cb, Cb> Cd(Bd.transpose() * Ad * Bd);

  auto baseVari
      = new quad_form_vari<std::remove_pointer_t<decltype(A_mem)>, Ra, Ca,
                           std::remove_pointer_t<decltype(B_mem)>, Rb, Cb>(
          Cd, A_mem, Ad_mem, A.rows(), B_mem, Bd_mem, B.rows(), B.cols(),
          false);

  return baseVari->getOutput();
}
}  // namespace

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
inline typename boost::enable_if_c<boost::is_same<Ta, var>::value
                                       || boost::is_same<Tb, var>::value,
                                   Eigen::Matrix<var, Cb, Cb> >::type
quad_form(const Eigen::Matrix<Ta, Ra, Ca>& A,
          const Eigen::Matrix<Tb, Rb, Cb>& B) {
  return quad_form_inner(A, B, false);
}

template <typename Ta, int Ra, int Ca, typename Tb, int Rb>
inline typename boost::enable_if_c<
    boost::is_same<Ta, var>::value || boost::is_same<Tb, var>::value, var>::type
quad_form(const Eigen::Matrix<Ta, Ra, Ca>& A,
          const Eigen::Matrix<Tb, Rb, 1>& B) {
  return quad_form_inner(A, B, false)(0, 0);
}

}  // namespace math
}  // namespace stan
#endif
