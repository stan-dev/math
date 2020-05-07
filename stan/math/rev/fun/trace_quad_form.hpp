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
template <typename EigMatL, typename EigMatR>
class trace_quad_form_vari_alloc : public chainable_alloc {
 public:
  trace_quad_form_vari_alloc(EigMatL&& A, EigMatR&& B)
      : A_(std::forward<EigMatL>(A)), B_(std::forward<EigMatR>(B)) {}

  double compute() { return trace_quad_form(value_of(A_), value_of(B_)); }

  typename std::decay_t<EigMatL>::PlainObject A_;
  typename std::decay_t<EigMatR>::PlainObject B_;
};

template <typename EigMatL, typename EigMatR>
class trace_quad_form_vari : public vari {
 protected:
  template <
      typename EigMatA, typename EigMatBd,
      require_all_eigen_vt<std::is_arithmetic, EigMatA, EigMatBd>* = nullptr>
  static inline void chainA(EigMatA&& A, EigMatBd&& Bd, double adjC) {}

  template <typename EigMatB, typename EigMatAd, typename EigMatBd,
            require_all_eigen_vt<std::is_arithmetic, EigMatB, EigMatAd,
                                 EigMatBd>* = nullptr>
  static inline void chainB(EigMatB&& B, EigMatAd&& Ad, EigMatBd&& Bd,
                            double adjC) {}

  template <typename EigMatA, typename EigMatBd,
            require_eigen_vt<is_var, EigMatA>* = nullptr,
            require_eigen_vt<std::is_arithmetic, EigMatBd>* = nullptr>
  static inline void chainA(EigMatA&& A, EigMatBd&& Bd, double adjC) {
    A.adj() += adjC * Bd * Bd.transpose();
  }

  template <
      typename EigMatB, typename EigMatAd, typename EigMatBd,
      require_eigen_vt<is_var, EigMatB>* = nullptr,
      require_all_eigen_vt<std::is_arithmetic, EigMatAd, EigMatBd>* = nullptr>
  static inline void chainB(EigMatB&& B, EigMatAd&& Ad, EigMatBd& Bd,
                            double adjC) {
    B.adj() += adjC * (Ad + Ad.transpose()) * Bd;
  }
  template <
      typename EigMatA, typename EigMatB, typename EigMatAd, typename EigMatBd,
      require_all_eigen_t<EigMatA, EigMatB>* = nullptr,
      require_all_eigen_vt<std::is_arithmetic, EigMatAd, EigMatBd>* = nullptr>
  inline void chainAB(EigMatA&& A, EigMatB&& B, EigMatAd&& Ad, EigMatBd&& Bd,
                      double adjC) {
    chainA(A, Bd, adjC);
    chainB(B, Ad, Bd, adjC);
  }

 public:
  explicit trace_quad_form_vari(
      trace_quad_form_vari_alloc<EigMatL, EigMatR>* impl)
      : vari(impl->compute()), impl_(impl) {}

  virtual void chain() {
    chainAB(impl_->A_, impl_->B_, value_of(impl_->A_), value_of(impl_->B_),
            adj_);
  }
  trace_quad_form_vari_alloc<EigMatL, EigMatR>* impl_;
};
}  // namespace internal

template <typename EigMatL, typename EigMatR,
          require_any_eigen_vt<is_var, EigMatL, EigMatR>* = nullptr>
inline auto trace_quad_form(EigMatL&& A, EigMatR&& B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);
  auto* baseVari = new internal::trace_quad_form_vari_alloc<EigMatL, EigMatR>(
      std::forward<EigMatL>(A), std::forward<EigMatR>(B));
  return var(new internal::trace_quad_form_vari<EigMatL, EigMatR>(baseVari));
}

}  // namespace math
}  // namespace stan
#endif
