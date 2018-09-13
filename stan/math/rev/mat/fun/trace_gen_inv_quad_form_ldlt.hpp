#ifndef STAN_MATH_REV_MAT_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_GEN_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <boost/utility/enable_if.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {

namespace {

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
class trace_gen_inv_quad_form_ldlt_vari : public vari {
  Td* D_mem_;
  const double* Dd_mem_;
  int D_rows_;
  int D_cols_;
  LDLT_factor<Ta, Ra, Ca> A_;
  Tb* B_mem_;
  int B_rows_;
  int B_cols_;
  const double* C_mem_;
  int C_rows_;
  int C_cols_;
  const double* AinvB_mem_;

 protected:
  static inline void chainD(const double* D_mem_,
                            Eigen::Map<const Eigen::Matrix<double, Cb, Cb> >& C,
                            double adj) {}

  static inline void chainA(
      const LDLT_factor<double, Ra, Ca>& A,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB,
      Eigen::Map<const Eigen::Matrix<double, Rd, Cd> >& Dd, double adj) {}

  static inline void chainB(
      const double* B_mem_,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB,
      Eigen::Map<const Eigen::Matrix<double, Rd, Cd> >& Dd, double adj) {}

  static inline void chainD(vari** D_mem_,
                            Eigen::Map<const Eigen::Matrix<double, Cb, Cb> >& C,
                            double adj) {
    for (int i = 0; i < C.size(); ++i)
      D_mem_[i]->adj_ += adj * C(i);
  }

  static inline void chainA(
      const LDLT_factor<var, Ra, Ca>& A,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB,
      Eigen::Map<const Eigen::Matrix<double, Rd, Cd> >& Dd, double adj) {
    Eigen::Matrix<double, Ra, Ca> adjA
        = -adj * (AinvB * Dd.transpose() * AinvB.transpose());

    for (int i = 0; i < adjA.size(); ++i)
      A.getVariA(i)->adj_ += adjA(i);
  }

  static inline void chainB(
      vari** B_mem_, Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB,
      Eigen::Map<const Eigen::Matrix<double, Rd, Cd> >& Dd, double adj) {
    Eigen::Matrix<double, Rb, Cb> adjB = adj * AinvB * (Dd + Dd.transpose());

    for (int i = 0; i < adjB.size(); i++)
      B_mem_[i]->adj_ += adjB(i);
  }

 public:
  explicit trace_gen_inv_quad_form_ldlt_vari(
      double value, Td* D_mem, const double* Dd_mem, int D_rows, int D_cols,
      const LDLT_factor<Ta, Ra, Ca>& A, Tb* B_mem, int B_rows, int B_cols,
      const double* C_mem, int C_rows, int C_cols, const double* AinvB_mem)
      : vari(value),
        D_mem_(D_mem),
        Dd_mem_(Dd_mem),
        D_rows_(D_rows),
        D_cols_(D_cols),
        A_(A),
        B_mem_(B_mem),
        B_rows_(B_rows),
        B_cols_(B_cols),
        C_mem_(C_mem),
        C_rows_(C_rows),
        C_cols_(C_cols),
        AinvB_mem_(AinvB_mem) {}

  virtual void chain() {
    // F = trace(D * B' * inv(A) * B)
    // aA = -aF * inv(A') * B * D' * B' * inv(A')
    // aB = aF*(inv(A) * B * D + inv(A') * B * D')
    // aD = aF*(B' * inv(A) * B)
    Eigen::Map<const Eigen::Matrix<double, Rd, Cd> > Dd(Dd_mem_, D_rows_,
                                                        D_cols_);
    Eigen::Map<const Eigen::Matrix<double, Rb, Cb> > AinvB(AinvB_mem_, B_rows_,
                                                           B_cols_);
    Eigen::Map<const Eigen::Matrix<double, Cb, Cb> > C(C_mem_, C_cols_,
                                                       C_cols_);
    chainD(D_mem_, C, adj_);
    chainA(A_, AinvB, Dd, adj_);
    chainB(B_mem_, AinvB, Dd, adj_);
  }
};

}  // namespace

/**
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       trace(D B^T A^-1 B)
 * where D is a square matrix and the LDLT_factor of A is provided.
 **/
template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
inline typename boost::enable_if_c<stan::is_var<Td>::value
                                       || stan::is_var<Ta>::value
                                       || stan::is_var<Tb>::value,
                                   var>::type
trace_gen_inv_quad_form_ldlt(const Eigen::Matrix<Td, Rd, Cd>& D,
                             const LDLT_factor<Ta, Ra, Ca>& A,
                             const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("trace_gen_inv_quad_form_ldlt", "D", D);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "A", A, "B", B);
  check_multiplicable("trace_gen_inv_quad_form_ldlt", "B", B, "D", D);

  auto D_mem = build_vari_pointer_array_if_necessary(D.data(), D.size());
  auto B_mem = build_vari_pointer_array_if_necessary(B.data(), B.size());

  const double* Dd_mem = build_double_array(D_mem, D.size());
  double* AinvB_mem = build_double_array(B_mem, B.size());
  Eigen::Map<Eigen::Matrix<double, Rb, Cb> > AinvB(AinvB_mem, B.rows(),
                                                   B.cols());
  Eigen::Matrix<double, Rb, Cb> Bd = AinvB;
  A.solveInPlace(AinvB);

  Eigen::Matrix<double, Ra, Ca> Dd
      = Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >(Dd_mem, D.rows(),
                                                         D.cols());

  double* C_mem = ChainableStack::instance().memalloc_.alloc_array<double>(
      B.cols() * B.cols());

  Eigen::Map<Eigen::Matrix<double, Cb, Cb> > C(C_mem, B.cols(), B.cols());
  C = Bd.transpose() * AinvB;

  double value = (Dd * C).trace();

  return var(new trace_gen_inv_quad_form_ldlt_vari<
             std::remove_pointer_t<decltype(D_mem)>, Rd, Cd, Ta, Ra, Ca,
             std::remove_pointer_t<decltype(B_mem)>, Rb, Cb>(
      value, D_mem, Dd_mem, D.rows(), D.cols(), A, B_mem, B.rows(), B.cols(),
      C_mem, C.rows(), C.cols(), AinvB_mem));
}

}  // namespace math
}  // namespace stan
#endif
