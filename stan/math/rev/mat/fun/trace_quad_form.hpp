#ifndef STAN_MATH_REV_MAT_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_QUAD_FORM_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/trace_quad_form.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {
namespace {
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_quad_form_vari : public vari {
 protected:
  Ta* A_mem_;
  const double* Ad_mem_;
  int A_rows_;
  Tb* B_mem_;
  const double* Bd_mem_;
  int B_rows_;
  int B_cols_;

  static inline void chainA(
      Eigen::Map<Eigen::Matrix<double, Ra, Ca> >& A,
      Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >& Ad, double adjC) {}
  static inline void chainB(
      Eigen::Map<Eigen::Matrix<double, Rb, Cb> >& B,
      Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >& Ad,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd, double adjC) {}

  static inline void chainA(
      Eigen::Map<Eigen::Matrix<var, Ra, Ca> >& A,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd, double adjC) {
    Eigen::Matrix<double, Ra, Ca> adjA(adjC * Bd * Bd.transpose());
    for (int i = 0; i < A.size(); ++i)
      A.data()[i].vi_->adj_ += adjA.data()[i];
  }
  static inline void chainB(
      Eigen::Map<Eigen::Matrix<var, Rb, Cb> >& B,
      Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >& Ad,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd, double adjC) {
    Eigen::Matrix<double, Ra, Ca> adjB(adjC * (Ad + Ad.transpose()) * Bd);
    for (int i = 0; i < B.size(); ++i)
      B.data()[i].vi_->adj_ += adjB.data()[i];
  }

 public:
  explicit trace_quad_form_vari(double result, Ta* A_mem, const double* Ad_mem,
                                int A_rows, Tb* B_mem, const double* Bd_mem,
                                int B_rows, int B_cols)
      : vari(result),
        A_mem_(A_mem),
        Ad_mem_(Ad_mem),
        A_rows_(A_rows),
        B_mem_(B_mem),
        Bd_mem_(Bd_mem),
        B_rows_(B_rows),
        B_cols_(B_cols) {}

  virtual void chain() {
    Eigen::Map<Eigen::Matrix<Ta, Ra, Ca> > A(A_mem_, A_rows_, A_rows_);
    Eigen::Map<Eigen::Matrix<Tb, Rb, Cb> > B(B_mem_, B_rows_, B_cols_);
    Eigen::Map<const Eigen::Matrix<double, Ra, Ca> > Ad(Ad_mem_, A_rows_,
                                                        A_rows_);
    Eigen::Map<const Eigen::Matrix<double, Rb, Cb> > Bd(Bd_mem_, B_rows_,
                                                        B_cols_);

    chainA(A, Bd, adj_);
    chainB(B, Ad, Bd, adj_);
  }
};
}  // namespace

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
inline typename boost::enable_if_c<
    boost::is_same<Ta, var>::value || boost::is_same<Tb, var>::value, var>::type
trace_quad_form(const Eigen::Matrix<Ta, Ra, Ca>& A,
                const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);

  int A_rows = A.rows();  // A is square
  int B_rows = B.rows();
  int B_cols = B.cols();
  Ta* A_mem = ChainableStack::instance().memalloc_.alloc_array<Ta>(A.size());
  Tb* B_mem = ChainableStack::instance().memalloc_.alloc_array<Tb>(B.size());

  for (int i = 0; i < A.size(); ++i)
    A_mem[i] = A.data()[i];
  for (int i = 0; i < B.size(); ++i)
    B_mem[i] = B.data()[i];

  double* Ad_mem = build_double_array_if_necessary(A_mem, A_rows * A_rows);
  double* Bd_mem = build_double_array_if_necessary(B_mem, B_rows * B_cols);

  Eigen::Matrix<double, Ra, Ca> Ad
      = Eigen::Map<Eigen::Matrix<double, Ra, Ca> >(Ad_mem, A_rows, A_rows);
  Eigen::Matrix<double, Ra, Ca> Bd
      = Eigen::Map<Eigen::Matrix<double, Rb, Cb> >(Bd_mem, B_rows, B_cols);

  double result = trace_quad_form(Ad, Bd);

  return var(new trace_quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>(
      result, A_mem, Ad_mem, A_rows, B_mem, Bd_mem, B_rows, B_cols));
}

}  // namespace math
}  // namespace stan
#endif
