#ifndef STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP
#define STAN_MATH_REV_MAT_FUN_TRACE_GEN_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/trace_gen_quad_form.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <type_traits>

namespace stan {
namespace math {
namespace internal {
template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
class trace_gen_quad_form_vari_alloc : public chainable_alloc {
 public:
  trace_gen_quad_form_vari_alloc(const Eigen::Matrix<Td, Rd, Cd>& D,
                                 const Eigen::Matrix<Ta, Ra, Ca>& A,
                                 const Eigen::Matrix<Tb, Rb, Cb>& B)
      : D_(D), A_(A), B_(B) {}

  double compute() {
    return trace_gen_quad_form(value_of(D_), value_of(A_), value_of(B_));
  }

  Eigen::Matrix<Td, Rd, Cd> D_;
  Eigen::Matrix<Ta, Ra, Ca> A_;
  Eigen::Matrix<Tb, Rb, Cb> B_;
};

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb>
class trace_gen_quad_form_vari : public vari {
 protected:
  static inline void computeAdjoints(double adj,
                                     const Eigen::Matrix<double, Rd, Cd>& D,
                                     const Eigen::Matrix<double, Ra, Ca>& A,
                                     const Eigen::Matrix<double, Rb, Cb>& B,
                                     Eigen::Matrix<var, Rd, Cd>* varD,
                                     Eigen::Matrix<var, Ra, Ca>* varA,
                                     Eigen::Matrix<var, Rb, Cb>* varB) {
    Eigen::Matrix<double, Ca, Cb> AtB;
    Eigen::Matrix<double, Ra, Cb> BD;
    if (varB || varA) {
      BD.noalias() = B * D;
    }
    if (varB || varD) {
      AtB.noalias() = A.transpose() * B;
    }

    if (varB) {
      (*varB).adj() += adj * (A * BD + AtB * D.transpose());
    }
    if (varA) {
      (*varA).adj() += adj * (B * BD.transpose());
    }
    if (varD) {
      (*varD).adj() += adj * (B.transpose() * AtB);
    }
  }

 public:
  explicit trace_gen_quad_form_vari(
      trace_gen_quad_form_vari_alloc<Td, Rd, Cd, Ta, Ra, Ca, Tb, Rb, Cb>* impl)
      : vari(impl->compute()), impl_(impl) {}

  virtual void chain() {
    computeAdjoints(adj_, value_of(impl_->D_), value_of(impl_->A_),
                    value_of(impl_->B_),
                    reinterpret_cast<Eigen::Matrix<var, Rd, Cd>*>(
                        std::is_same<Td, var>::value ? (&impl_->D_) : NULL),
                    reinterpret_cast<Eigen::Matrix<var, Ra, Ca>*>(
                        std::is_same<Ta, var>::value ? (&impl_->A_) : NULL),
                    reinterpret_cast<Eigen::Matrix<var, Rb, Cb>*>(
                        std::is_same<Tb, var>::value ? (&impl_->B_) : NULL));
  }

  trace_gen_quad_form_vari_alloc<Td, Rd, Cd, Ta, Ra, Ca, Tb, Rb, Cb>* impl_;
};
}  // namespace internal

template <typename Td, int Rd, int Cd, typename Ta, int Ra, int Ca, typename Tb,
          int Rb, int Cb, typename = require_any_var_t<Td, Ta, Tb>>
inline var trace_gen_quad_form(const Eigen::Matrix<Td, Rd, Cd>& D,
                               const Eigen::Matrix<Ta, Ra, Ca>& A,
                               const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("trace_gen_quad_form", "A", A);
  check_square("trace_gen_quad_form", "D", D);
  check_multiplicable("trace_gen_quad_form", "A", A, "B", B);
  check_multiplicable("trace_gen_quad_form", "B", B, "D", D);

  internal::trace_gen_quad_form_vari_alloc<Td, Rd, Cd, Ta, Ra, Ca, Tb, Rb, Cb>*
      baseVari
      = new internal::trace_gen_quad_form_vari_alloc<Td, Rd, Cd, Ta, Ra, Ca, Tb,
                                                     Rb, Cb>(D, A, B);

  return var(new internal::trace_gen_quad_form_vari<Td, Rd, Cd, Ta, Ra, Ca, Tb,
                                                    Rb, Cb>(baseVari));
}

}  // namespace math
}  // namespace stan
#endif
