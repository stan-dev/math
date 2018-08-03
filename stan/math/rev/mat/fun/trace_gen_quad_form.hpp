#ifndef STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/trace_gen_quad_form.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {
namespace {
template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
class trace_gen_quad_form_vari : public vari {
 protected:
  Td* D_mem_;
  const double* Dd_mem_;
  int D_rows_;
  Ta* A_mem_;
  const double* Ad_mem_;
  int A_rows_;
  Tb* B_mem_;
  const double* Bd_mem_;
  int B_rows_;
  int B_cols_;

  static inline void computeAdjoints(
      double adj, Eigen::Map<const Eigen::Matrix<double, Cd, Rd> >& D,
      Eigen::Map<const Eigen::Matrix<double, Ca, Ra> >& A,
      Eigen::Map<const Eigen::Matrix<double, Cb, Rb> >& B, vari** D_mem,
      vari** A_mem, vari** B_mem) {
    Eigen::Matrix<double, Ca, Cb> AtB;
    Eigen::Matrix<double, Ra, Cb> BD;

    if (B_mem || A_mem)
      BD.noalias() = B * D;
    if (B_mem || D_mem)
      AtB.noalias() = A.transpose() * B;

    if (D_mem) {
      Eigen::Matrix<double, Rd, Cd> adjD(adj * (B.transpose() * AtB));
      for (int i = 0; i < adjD.size(); ++i)
        D_mem[i]->adj_ += adjD.data()[i];
    }

    if (A_mem) {
      Eigen::Matrix<double, Ra, Ca> adjA(adj * (B * BD.transpose()));
      for (int i = 0; i < adjA.size(); ++i)
        A_mem[i]->adj_ += adjA.data()[i];
    }

    if (B_mem) {
      Eigen::Matrix<double, Rb, Cb> adjB(adj * (A * BD + AtB * D.transpose()));
      for (int i = 0; i < adjB.size(); ++i)
        B_mem[i]->adj_ += adjB.data()[i];
    }
  }

  vari** null_if_not_vari_pointer_array(vari** ptr) { return ptr; }

  vari** null_if_not_vari_pointer_array(const double* ptr) {
    return static_cast<vari**>(NULL);
  }

 public:
  explicit trace_gen_quad_form_vari(double result, Td* D_mem,
                                    const double* Dd_mem, int D_rows, Ta* A_mem,
                                    const double* Ad_mem, int A_rows, Tb* B_mem,
                                    const double* Bd_mem, int B_rows,
                                    int B_cols)
      : vari(result),
        D_mem_(D_mem),
        Dd_mem_(Dd_mem),
        D_rows_(D_rows),
        A_mem_(A_mem),
        Ad_mem_(Ad_mem),
        A_rows_(A_rows),
        B_mem_(B_mem),
        Bd_mem_(Bd_mem),
        B_rows_(B_rows),
        B_cols_(B_cols) {}

  virtual void chain() {
    Eigen::Map<const Eigen::Matrix<double, Cd, Rd> > D(Dd_mem_, D_rows_,
                                                       D_rows_);
    Eigen::Map<const Eigen::Matrix<double, Ca, Ra> > A(Ad_mem_, A_rows_,
                                                       A_rows_);
    Eigen::Map<const Eigen::Matrix<double, Cb, Rb> > B(Bd_mem_, B_rows_,
                                                       B_cols_);

    computeAdjoints(adj_, D, A, B, null_if_not_vari_pointer_array(D_mem_),
                    null_if_not_vari_pointer_array(A_mem_),
                    null_if_not_vari_pointer_array(B_mem_));
  }
};
}  // namespace

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
inline typename boost::enable_if_c<boost::is_same<Td, var>::value
                                       || boost::is_same<Ta, var>::value
                                       || boost::is_same<Tb, var>::value,
                                   var>::type
trace_gen_quad_form(const Eigen::Matrix<Td, Rd, Cd>& D,
                    const Eigen::Matrix<Ta, Ra, Ca>& A,
                    const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);

  auto D_mem = build_vari_pointer_array_if_necessary(D.data(), D.size());
  auto A_mem = build_vari_pointer_array_if_necessary(A.data(), A.size());
  auto B_mem = build_vari_pointer_array_if_necessary(B.data(), B.size());

  const double* Dd_mem = build_double_array(D_mem, D.size());
  const double* Ad_mem = build_double_array(A_mem, A.size());
  const double* Bd_mem = build_double_array(B_mem, B.size());

  Eigen::Matrix<double, Rd, Cd> Dd
      = Eigen::Map<const Eigen::Matrix<double, Rd, Cd> >(Dd_mem, D.rows(),
                                                         D.cols());
  Eigen::Matrix<double, Ra, Ca> Ad
      = Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >(Ad_mem, A.rows(),
                                                         A.cols());
  Eigen::Matrix<double, Rb, Cb> Bd
      = Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >(Bd_mem, B.rows(),
                                                         B.cols());

  double result = trace_gen_quad_form(Dd, Ad, Bd);

  auto baseVari = new trace_gen_quad_form_vari<
      std::remove_pointer_t<decltype(D_mem)>, Rd, Cd,
      std::remove_pointer_t<decltype(A_mem)>, Ra, Ca,
      std::remove_pointer_t<decltype(B_mem)>, Rb, Cb>(
      result, D_mem, Dd_mem, D.rows(), A_mem, Ad_mem, A.rows(), B_mem, Bd_mem,
      B.rows(), B.cols());

  return var(baseVari);
}

}  // namespace math
}  // namespace stan
#endif
