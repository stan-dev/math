#ifndef STAN_MATH_REV_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_TRACE_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/trace_quad_form.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <type_traits>

namespace stan {
namespace math {
namespace internal {
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_quad_form_vari_alloc : public chainable_alloc {
 public:
  trace_quad_form_vari_alloc(const Eigen::Matrix<Ta, Ra, Ca>& A,
                             const Eigen::Matrix<Tb, Rb, Cb>& B)
      : A_(A), B_(B) {}

  double compute() { return trace_quad_form(value_of(A_), value_of(B_)); }

  Eigen::Matrix<Ta, Ra, Ca> A_;
  Eigen::Matrix<Tb, Rb, Cb> B_;
};

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class trace_quad_form_vari : public vari {
 protected:
  static inline void chainA(Eigen::Matrix<double, Ra, Ca>& A,
                            const Eigen::Matrix<double, Rb, Cb>& Bd,
                            double adjC) {}
  static inline void chainB(Eigen::Matrix<double, Rb, Cb>& B,
                            const Eigen::Matrix<double, Ra, Ca>& Ad,
                            const Eigen::Matrix<double, Rb, Cb>& Bd,
                            double adjC) {}

  static inline void chainA(Eigen::Matrix<var, Ra, Ca>& A,
                            const Eigen::Matrix<double, Rb, Cb>& Bd,
                            double adjC) {
    A.adj() += adjC * Bd * Bd.transpose();
  }
  static inline void chainB(Eigen::Matrix<var, Rb, Cb>& B,
                            const Eigen::Matrix<double, Ra, Ca>& Ad,
                            const Eigen::Matrix<double, Rb, Cb>& Bd,
                            double adjC) {
    B.adj() += adjC * (Ad + Ad.transpose()) * Bd;
  }

  inline void chainAB(Eigen::Matrix<Ta, Ra, Ca>& A,
                      Eigen::Matrix<Tb, Rb, Cb>& B,
                      const Eigen::Matrix<double, Ra, Ca>& Ad,
                      const Eigen::Matrix<double, Rb, Cb>& Bd, double adjC) {
    chainA(A, Bd, adjC);
    chainB(B, Ad, Bd, adjC);
  }

 public:
  explicit trace_quad_form_vari(
      trace_quad_form_vari_alloc<Ta, Ra, Ca, Tb, Rb, Cb>* impl)
      : vari(impl->compute()), impl_(impl) {}

  virtual void chain() {
    chainAB(impl_->A_, impl_->B_, value_of(impl_->A_), value_of(impl_->B_),
            adj_);
  }

  trace_quad_form_vari_alloc<Ta, Ra, Ca, Tb, Rb, Cb>* impl_;
};
}  // namespace internal

template <typename EigMat1, typename EigMat2,
          typename = require_any_vt_var<EigMat1, EigMat2>>
inline return_type_t<EigMat1, EigMat2> trace_quad_form(const EigMat1& A,
                                                       const EigMat2& B) {
  using Ta = value_type_t<EigMat1>;
  using Tb = value_type_t<EigMat2>;
  constexpr int Ra = EigMat1::RowsAtCompileTime;
  constexpr int Ca = EigMat1::ColsAtCompileTime;
  constexpr int Rb = EigMat2::RowsAtCompileTime;
  constexpr int Cb = EigMat2::ColsAtCompileTime;
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);

  auto* baseVari
      = new internal::trace_quad_form_vari_alloc<Ta, Ra, Ca, Tb, Rb, Cb>(A, B);

  return var(
      new internal::trace_quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>(baseVari));
}

}  // namespace math
}  // namespace stan
#endif
