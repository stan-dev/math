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
      const double* A_mem, int A_rows_,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd, double adjC) {}
  static inline void chainB(
      const double* B_mem, Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >& Ad,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd, double adjC) {}

  static inline void chainA(
      vari** A_mem, int A_rows_,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd, double adjC) {
    if (Bd.cols()
        == 1) {  // We can easily avoid building a full matrix if B is a vector
      for (int j = 0; j < A_rows_; ++j)
        for (int i = 0; i < A_rows_; ++i)
          A_mem[i + j * A_rows_]->adj_ += adjC * Bd(i) * Bd(j);
    } else {
      Eigen::Matrix<double, Rb, Rb> adjB(adjC * Bd * Bd.transpose());
      for (int i = 0; i < adjB.size(); ++i)
        A_mem[i]->adj_ += adjB(i);
    }
  }
  static inline void chainB(
      vari** B_mem, Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >& Ad,
      Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >& Bd, double adjC) {
    Eigen::Matrix<double, Ra, Ca> adjB(adjC * (Ad + Ad.transpose()) * Bd);
    for (int i = 0; i < adjB.size(); ++i)
      B_mem[i]->adj_ += adjB.data()[i];
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
    Eigen::Map<const Eigen::Matrix<double, Ra, Ca> > Ad(Ad_mem_, A_rows_,
                                                        A_rows_);
    Eigen::Map<const Eigen::Matrix<double, Rb, Cb> > Bd(Bd_mem_, B_rows_,
                                                        B_cols_);

    chainA(A_mem_, A_rows_, Bd, adj_);
    chainB(B_mem_, Ad, Bd, adj_);
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

  auto A_mem = build_vari_pointer_array_if_necessary(A.data(), A.size());
  auto B_mem = build_vari_pointer_array_if_necessary(B.data(), B.size());

  const double* Ad_mem = build_double_array(A_mem, A.size());
  const double* Bd_mem = build_double_array(B_mem, B.size());

  Eigen::Matrix<double, Ra, Ca> Ad
      = Eigen::Map<const Eigen::Matrix<double, Ra, Ca> >(Ad_mem, A.rows(),
                                                         A.cols());
  Eigen::Matrix<double, Rb, Cb> Bd
      = Eigen::Map<const Eigen::Matrix<double, Rb, Cb> >(Bd_mem, B.rows(),
                                                         B.cols());

  double result = trace_quad_form(Ad, Bd);

  auto baseVari
      = new trace_quad_form_vari<std::remove_pointer_t<decltype(A_mem)>, Ra, Ca,
                                 std::remove_pointer_t<decltype(B_mem)>, Rb,
                                 Cb>(result, A_mem, Ad_mem, A.rows(), B_mem,
                                     Bd_mem, B.rows(), B.cols());

  return var(baseVari);
}

}  // namespace math
}  // namespace stan
#endif
