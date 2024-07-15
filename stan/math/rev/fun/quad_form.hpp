#ifndef STAN_MATH_REV_FUN_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_var_value.hpp>
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
    if (sym_) {
      matrix_d M = 0.5 * (Cd + Cd.transpose());
      Cd = M;
    }
    for (int j = 0; j < C_.cols(); j++) {
      for (int i = 0; i < C_.rows(); i++) {
        C_(i, j) = var(new vari(Cd(i, j), false));
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

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam Mat1 type of the first (square) matrix
 * @tparam Mat2 type of the second matrix
 *
 * @param A square matrix
 * @param B second matrix
 * @param symmetric indicates whether the output should be made symmetric
 * @return The quadratic form
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_var_matrix_t<Mat1, Mat2>* = nullptr>
inline auto quad_form_impl(Mat1&& A, Mat2&& B, bool symmetric) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  using return_t
      = return_var_matrix_t<decltype(value_of(B).transpose().eval()
                                     * value_of(A) * value_of(B).eval()),
                            Mat1, Mat2>;

  arena_t<Mat1> arena_A = std::forward<Mat1>(A);
  arena_t<Mat2> arena_B = std::forward<Mat2>(B);
  if constexpr (is_autodiffable_v<Mat1, Mat2>) {
    check_not_nan("multiply", "A", arena_A.val());
    check_not_nan("multiply", "B", arena_B.val());

    auto res_vals = to_arena(arena_B.val_op().transpose() * arena_A.val_op()
                             * arena_B.val_op());

    if (symmetric) {
      res_vals += res_vals.transpose().eval();
    }
    arena_t<return_t> res = std::move(res_vals);
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      auto C_adj_B_t = (res.adj() * arena_B.val_op().transpose()).eval();

      if constexpr (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += arena_B.val_op() * C_adj_B_t;
      } else {
        arena_A.adj() += arena_B.val_op() * C_adj_B_t;
      }

      if constexpr (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias()
            += arena_A.val_op() * C_adj_B_t.transpose()
               + arena_A.val_op().transpose() * arena_B.val_op() * res.adj();
      } else {
        arena_B.adj()
            += arena_A.val_op() * C_adj_B_t.transpose()
               + arena_A.val_op().transpose() * arena_B.val_op() * res.adj();
      }
    });

    return res;
  } else if constexpr (is_autodiffable_v<Mat2>) {
    check_not_nan("multiply", "A", arena_A);
    check_not_nan("multiply", "B", arena_B.val());

    auto res_vals
        = to_arena(arena_B.val_op().transpose() * arena_A * arena_B.val_op());

    if (symmetric) {
      res_vals += res_vals.transpose().eval();
    }
    arena_t<return_t> res = std::move(res_vals);
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      auto C_adj_B_t = (res.adj() * arena_B.val_op().transpose());

      if (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias()
            += arena_A * C_adj_B_t.transpose()
               + arena_A.transpose() * arena_B.val_op() * res.adj();
      } else {
        arena_B.adj() += arena_A * C_adj_B_t.transpose()
                         + arena_A.transpose() * arena_B.val_op() * res.adj();
      }
    });

    return res;
  } else if constexpr (is_autodiffable_v<Mat1>) {
    check_not_nan("multiply", "A", arena_A.val());
    check_not_nan("multiply", "B", arena_B);

    auto res_vals = to_arena(arena_B.transpose() * arena_A.val() * arena_B);

    if (symmetric) {
      res_vals += res_vals.transpose().eval();
    }
    arena_t<return_t> res = std::move(res_vals);
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      auto C_adj_B_t = (res.adj() * arena_B.transpose());

      if constexpr (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += arena_B * C_adj_B_t;
      } else {
        arena_A.adj() += arena_B * C_adj_B_t;
      }
    });

    return res;
  }
}
}  // namespace internal

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam EigMat1 type of the first (square) matrix
 * @tparam EigMat2 type of the second matrix
 *
 * @param A square matrix
 * @param B second matrix
 * @param symmetric indicates whether the output should be made symmetric
 * @return The quadratic form, which is a symmetric matrix.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename EigMat1, typename EigMat2,
          require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
          require_not_eigen_col_vector_t<EigMat2>* = nullptr,
          require_any_vt_var<EigMat1, EigMat2>* = nullptr>
inline promote_scalar_t<var, EigMat2> quad_form(const EigMat1& A,
                                                const EigMat2& B,
                                                bool symmetric = false) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  auto* baseVari = new internal::quad_form_vari<
      value_type_t<EigMat1>, EigMat1::RowsAtCompileTime,
      EigMat1::ColsAtCompileTime, value_type_t<EigMat2>,
      EigMat2::RowsAtCompileTime, EigMat2::ColsAtCompileTime>(A, B, symmetric);

  return baseVari->impl_->C_;
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * @tparam EigMat type of the matrix
 * @tparam ColVec type of the vector
 *
 * @param A square matrix
 * @param B vector
 * @param symmetric indicates whether the output should be made symmetric
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename EigMat, typename ColVec, require_eigen_t<EigMat>* = nullptr,
          require_eigen_col_vector_t<ColVec>* = nullptr,
          require_any_vt_var<EigMat, ColVec>* = nullptr>
inline var quad_form(const EigMat& A, const ColVec& B, bool symmetric = false) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  auto* baseVari = new internal::quad_form_vari<
      value_type_t<EigMat>, EigMat::RowsAtCompileTime,
      EigMat::ColsAtCompileTime, value_type_t<ColVec>,
      ColVec::RowsAtCompileTime, 1>(A, B, symmetric);

  return baseVari->impl_->C_(0, 0);
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * This overload handles arguments where one of Mat1 or Mat2 are
 * `var_value<T>` where `T` is an Eigen type. The other type can
 * also be a `var_value` or it can be a matrix type that inherits
 * from EigenBase
 *
 * @tparam Mat1 type of the first (square) matrix
 * @tparam Mat2 type of the second matrix
 *
 * @param A square matrix
 * @param B second matrix
 * @param symmetric indicates whether the output should be made symmetric
 * @return The quadratic form, which is a symmetric matrix.
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_not_col_vector_t<Mat2>* = nullptr,
          require_any_var_matrix_t<Mat1, Mat2>* = nullptr>
inline auto quad_form(Mat1&& A, Mat2&& B, bool symmetric = false) {
  return internal::quad_form_impl(std::forward<Mat1>(A), std::forward<Mat2>(B),
                                  symmetric);
}

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * This overload handles arguments where one of Mat or Vec are
 * `var_value<T>` where `T` is an Eigen type. The other type can
 * also be a `var_value`, or it can be a type that inherits from
 * EigenBase
 *
 * @tparam Mat type of the matrix
 * @tparam Vec type of the vector
 *
 * @param A square matrix
 * @param B vector
 * @param symmetric indicates whether the output should be made symmetric
 * @return The quadratic form (a scalar).
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Mat, typename Vec, require_matrix_t<Mat>* = nullptr,
          require_col_vector_t<Vec>* = nullptr,
          require_any_var_matrix_t<Mat, Vec>* = nullptr>
inline var quad_form(Mat&& A, Vec&& B, bool symmetric = false) {
  return internal::quad_form_impl(std::forward<Mat>(A), std::forward<Vec>(B),
                                  symmetric)(0, 0);
}

}  // namespace math
}  // namespace stan
#endif
