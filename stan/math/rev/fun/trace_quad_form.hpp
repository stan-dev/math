#ifndef STAN_MATH_REV_FUN_TRACE_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_TRACE_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/fun/to_var_value.hpp>
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
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_any_st_var<EigMat1, EigMat2>* = nullptr>
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

/**
 * Compute trace(B^T A B).
 *
 * This overload handles arguments where one of Mat1 or Mat2 are
 * `var_value<T>` where `T` is an Eigen type. The other type can
 * also be a `var_value` or it can be a type that inherits
 * from EigenBase
 *
 * @tparam Mat1 type of the first matrix
 * @tparam Mat2 type of the second matrix
 *
 * @param A matrix
 * @param B matrix
 * @return The trace of B^T A B
 * @throw std::domain_error if A is not square
 * @throw std::domain_error if A cannot be multiplied by B
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_var_matrix_t<Mat1, Mat2>* = nullptr>
inline var trace_quad_form(const Mat1& A, const Mat2& B) {
  check_square("trace_quad_form", "A", A);
  check_multiplicable("trace_quad_form", "A", A, "B", B);

  var res;

  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_A = A;
    arena_t<promote_scalar_t<var, Mat2>> arena_B = B;

    res = (value_of(arena_B).transpose() * value_of(arena_A)
           * value_of(arena_B))
              .trace();

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias()
            += res.adj() * value_of(arena_B) * value_of(arena_B).transpose();
      } else {
        arena_A.adj()
            += res.adj() * value_of(arena_B) * value_of(arena_B).transpose();
      }

      if (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias()
            += res.adj() * (value_of(arena_A) + value_of(arena_A).transpose())
               * value_of(arena_B);
      } else {
        arena_B.adj() += res.adj()
                         * (value_of(arena_A) + value_of(arena_A).transpose())
                         * value_of(arena_B);
      }
    });
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_A = value_of(A);
    arena_t<promote_scalar_t<var, Mat2>> arena_B = B;

    res = (value_of(arena_B).transpose() * value_of(arena_A)
           * value_of(arena_B))
              .trace();

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias()
            += res.adj() * (arena_A + arena_A.transpose()) * value_of(arena_B);
      } else {
        arena_B.adj()
            += res.adj() * (arena_A + arena_A.transpose()) * value_of(arena_B);
      }
    });
  } else {
    arena_t<promote_scalar_t<var, Mat1>> arena_A = A;
    arena_t<promote_scalar_t<double, Mat2>> arena_B = value_of(B);

    res = (arena_B.transpose() * value_of(arena_A) * arena_B).trace();

    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += res.adj() * arena_B * arena_B.transpose();
      } else {
        arena_A.adj() += res.adj() * arena_B * arena_B.transpose();
      }
    });
  }

  return res;
}

}  // namespace math
}  // namespace stan
#endif
