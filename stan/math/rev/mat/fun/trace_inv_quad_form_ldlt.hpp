#ifndef STAN_MATH_REV_MAT_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_INV_QUAD_FORM_LDLT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/LDLT_factor.hpp>
#include <boost/utility/enable_if.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>

namespace stan {
namespace math {

namespace {
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_inv_quad_form_ldlt_vari : public vari {
  LDLT_factor<Ta, Ra, Ca> A_;
  Tb* B_mem_;
  int B_rows_;
  int B_cols_;
  const double* AinvB_mem_;

 protected:
  static inline void chainA(
      const LDLT_factor<double, Ra, Ca>& A,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB, double adj) {}

  static inline void chainB(
      const double* B_mem_,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB, double adj) {}

  static inline void chainA(
      const LDLT_factor<var, Ra, Ca>& A,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB, double adj) {
    Eigen::Matrix<double, Ra, Ca> adjA = -adj * (AinvB * AinvB.transpose());

    for (int i = 0; i < adjA.size(); ++i)
      A.getVariA(i)->adj_ += adjA(i);
  }

  static inline void chainB(
      vari** B_mem_, Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& AinvB,
      double adj) {
    Eigen::Matrix<double, Rb, Cb> adjB = 2 * adj * AinvB;

    for (int i = 0; i < adjB.size(); i++)
      B_mem_[i]->adj_ += adjB(i);
  }

 public:
  explicit trace_inv_quad_form_ldlt_vari(double value,
                                         const LDLT_factor<Ta, Ra, Ca>& A,
                                         Tb* B_mem, int B_rows, int B_cols,
                                         const double* AinvB_mem)
      : vari(value),
        A_(A),
        B_mem_(B_mem),
        B_rows_(B_rows),
        B_cols_(B_cols),
        AinvB_mem_(AinvB_mem) {}

  virtual void chain() {
    // F = trace(B' * inv(A) * B)
    // aA = -aF * inv(A') * B * B' * inv(A')
    // aB = aF*(2 * inv(A) * B)
    Eigen::Map<const Eigen::Matrix<double, Rb, Cb> > AinvB(AinvB_mem_, B_rows_,
                                                           B_cols_);
    chainA(A_, AinvB, adj_);
    chainB(B_mem_, AinvB, adj_);
  }
};

}  // namespace

/**
 * Compute the trace of an inverse quadratic form.  I.E., this computes
 *       trace(B^T A^-1 B)
 * where the LDLT_factor of A is provided.
 **/
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
inline typename boost::enable_if_c<
    stan::is_var<Ta>::value || stan::is_var<Tb>::value, var>::type
trace_inv_quad_form_ldlt(const LDLT_factor<Ta, Ra, Ca>& A,
                         const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_multiplicable("trace_inv_quad_form_ldlt", "A", A, "B", B);

  auto B_mem = build_vari_pointer_array_if_necessary(B.data(), B.size());

  double* AinvB_mem = build_double_array(B_mem, B.size());
  Eigen::Map<Eigen::Matrix<double, Rb, Cb> > AinvB(AinvB_mem, B.rows(),
                                                   B.cols());
  Eigen::Matrix<double, Rb, Cb> Bd = AinvB;
  A.solveInPlace(AinvB);

  double value = (Bd.transpose() * AinvB).trace();

  return var(new trace_inv_quad_form_ldlt_vari<
             Ta, Ra, Ca, std::remove_pointer_t<decltype(B_mem)>, Rb, Cb>(
      value, A, B_mem, B.rows(), B.cols(), AinvB_mem));
}

}  // namespace math
}  // namespace stan
#endif
