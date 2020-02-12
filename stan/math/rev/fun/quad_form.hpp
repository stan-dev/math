#ifndef STAN_MATH_REV_FUN_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/quad_form.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class quad_form_vari_alloc : public chainable_alloc {
 private:
  inline void compute(const Eigen::Matrix<double, Ra, Ca>& A,
                      const Eigen::Matrix<double, Rb, Cb>& B) {
    matrix_d Cd = B.transpose() * A * B;
    for (int j = 0; j < C_.cols(); j++) {
      for (int i = 0; i < C_.rows(); i++) {
        if (sym_) {
          C_(i, j) = var(new vari(0.5 * (Cd(i, j) + Cd(j, i)), false));
        } else {
          C_(i, j) = var(new vari(Cd(i, j), false));
        }
      }
    }
  }

 public:
  quad_form_vari_alloc(const Eigen::Matrix<Ta, Ra, Ca>& A,
                       const Eigen::Matrix<Tb, Rb, Cb>& B,
                       bool symmetric = false)
      : A_(A), B_(B), C_(B_.cols(), B_.cols()), sym_(symmetric) {
    compute(value_of(A), value_of(B));
  }

  Eigen::Matrix<Ta, Ra, Ca> A_;
  Eigen::Matrix<Tb, Rb, Cb> B_;
  Eigen::Matrix<var, Cb, Cb> C_;
  bool sym_;
};

template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb>
class quad_form_vari : public vari {
 protected:
  inline void chainA(Eigen::Matrix<double, Ra, Ca>& A,
                     const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {}
  inline void chainB(Eigen::Matrix<double, Rb, Cb>& B,
                     const Eigen::Matrix<double, Ra, Ca>& Ad,
                     const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {}

  inline void chainA(Eigen::Matrix<var, Ra, Ca>& A,
                     const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {
    A.adj() += Bd * adjC * Bd.transpose();
  }
  inline void chainB(Eigen::Matrix<var, Rb, Cb>& B,
                     const Eigen::Matrix<double, Ra, Ca>& Ad,
                     const Eigen::Matrix<double, Rb, Cb>& Bd,
                     const Eigen::Matrix<double, Cb, Cb>& adjC) {
    B.adj() += Ad * Bd * adjC.transpose() + Ad.transpose() * Bd * adjC;
  }

  inline void chainAB(Eigen::Matrix<Ta, Ra, Ca>& A,
                      Eigen::Matrix<Tb, Rb, Cb>& B,
                      const Eigen::Matrix<double, Ra, Ca>& Ad,
                      const Eigen::Matrix<double, Rb, Cb>& Bd,
                      const Eigen::Matrix<double, Cb, Cb>& adjC) {
    chainA(A, Bd, adjC);
    chainB(B, Ad, Bd, adjC);
  }

 public:
  quad_form_vari(const Eigen::Matrix<Ta, Ra, Ca>& A,
                 const Eigen::Matrix<Tb, Rb, Cb>& B, bool symmetric = false)
      : vari(0.0) {
    impl_ = new quad_form_vari_alloc<Ta, Ra, Ca, Tb, Rb, Cb>(A, B, symmetric);
  }

  virtual void chain() {
    matrix_d adjC = impl_->C_.adj();

    chainAB(impl_->A_, impl_->B_, value_of(impl_->A_), value_of(impl_->B_),
            adjC);
  }

  quad_form_vari_alloc<Ta, Ra, Ca, Tb, Rb, Cb>* impl_;
};
}  // namespace internal

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam Ta type of elements in the square matrix
 * @tparam Ra number of rows in the square matrix, can be Eigen::Dynamic
 * @tparam Ca number of columns in the square matrix, can be Eigen::Dynamic
 * @tparam Tb type of elements in the second matrix
 * @tparam Rb number of rows in the second matrix, can be Eigen::Dynamic
 * @tparam Cb number of columns in the second matrix, can be Eigen::Dynamic
 *
 * @param A square matrix
 * @param B second matrix
 * @return The quadratic form, which is a symmetric matrix of size Cb.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Rb, int Cb,
          require_any_var_t<Ta, Tb>...>
inline Eigen::Matrix<var, Cb, Cb> quad_form(
    const Eigen::Matrix<Ta, Ra, Ca>& A, const Eigen::Matrix<Tb, Rb, Cb>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>* baseVari
      = new internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, Cb>(A, B);

  return baseVari->impl_->C_;
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * @tparam Ta type of elements in the square matrix
 * @tparam Ra number of rows in the square matrix, can be Eigen::Dynamic
 * @tparam Ca number of columns in the square matrix, can be Eigen::Dynamic
 * @tparam Tb type of elements in the vector
 * @tparam Rb number of rows in the vector, can be Eigen::Dynamic
 *
 * @param A square matrix
 * @param B vector
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Ta, int Ra, int Ca, typename Tb, int Rb,
          require_any_var_t<Ta, Tb>...>
inline var quad_form(const Eigen::Matrix<Ta, Ra, Ca>& A,
                     const Eigen::Matrix<Tb, Rb, 1>& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>* baseVari
      = new internal::quad_form_vari<Ta, Ra, Ca, Tb, Rb, 1>(A, B);

  return baseVari->impl_->C_(0, 0);
}

}  // namespace math
}  // namespace stan
#endif
