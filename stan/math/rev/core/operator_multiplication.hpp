#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/op_vari.hpp>
#include <stan/math/rev/functor/adj_jac_apply.hpp>
#include <stan/math/rev/meta/is_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isnan.hpp>

namespace stan {
namespace math {

namespace internal {
/**
 * Base for multiplication, to be specliazed for chain types.
 */
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class multiply_vari {};

/**
 * Specialization of var multiplication for two `var_value`.
 */
template <typename VariVal, typename Vari1, typename Vari2>
class multiply_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>>
    final : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;

 public:
  multiply_vari(Vari1* avi, Vari2* bvi)
      : op_vari<VariVal, Vari1*, Vari2*>(avi->val_ * bvi->val_, avi, bvi) {}
  /**
   * `chain_impl` is called from `chain()` and exists so one specialized struct
   * can call either the scalar or matrix `chain()` methods. SFINAE only works
   * on "deduced" template types. So the trick here is to make template types,
   * T1 and T2, set their defaults to the class template types, then do the
   * regular `requires`. Since `chain_impl` has no inputs to deduce the
   * template types will always fall back to their default values. Since the
   * compiler has "deduced" these types we can these use the standard requires
   * to SFINAE out either the arithmetic or matrix version.
   */
  template <typename T1 = Vari1, typename T2 = Vari2,
            require_all_vari_vt<std::is_arithmetic, T1, T2>* = nullptr>
  inline void chain_impl() {
    avi()->adj_ += bvi()->val_ * this->adj_;
    bvi()->adj_ += avi()->val_ * this->adj_;
  }

  template <typename T1 = Vari1, typename T2 = Vari2,
            require_all_vari_vt<is_eigen, T1, T2>* = nullptr>
  inline void chain_impl() {
    avi()->adj_ += this->adj_ * bvi()->val_.transpose();
    bvi()->adj_ += avi()->val_.transpose() * this->adj_;
  }

  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      fill(avi()->adj_, NOT_A_NUMBER);
      fill(bvi()->adj_, NOT_A_NUMBER);
    } else {
      chain_impl();
    }
  }
};

/**
 * Specialization of var multiplication for `var_value` and arithmetic
 */
template <typename VariVal, typename Vari, typename Arith>
class multiply_vari<VariVal, Vari, Arith, require_vt_arithmetic<Arith>> final
    : public op_vari<VariVal, Vari*, Arith> {
  using op_vari<VariVal, Vari*, Arith>::avi;
  using op_vari<VariVal, Vari*, Arith>::bd;

 public:
  multiply_vari(Vari* avi, const Arith& b)
      : op_vari<VariVal, Vari*, Arith>(avi->val_ * b, avi, b) {}

  template <typename T1 = Vari, typename T2 = Arith,
            require_vari_vt<std::is_arithmetic, T1>* = nullptr,
            require_arithmetic_t<T2>* = nullptr>
  inline void chain_impl() {
    avi()->adj_ += this->adj_ * bd();
  }

  template <typename T1 = Vari, typename T2 = Arith,
            require_vari_vt<is_eigen, T1>* = nullptr,
            require_vt_arithmetic<T2>* = nullptr>
  inline void chain_impl() {
    avi()->adj_ += this->adj_ * bd().transpose();
  }
  // NOTE: THIS IS WRONG
  template <typename T1 = Arith, typename T2 = Vari,
            require_eigen_t<T1>* = nullptr,
            require_vari_vt<std::is_arithmetic, T2>* = nullptr>
  void chain_impl() {
    avi()->adj_ += this->adj_.sum();
  }

  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      fill(avi()->adj_, NOT_A_NUMBER);
    } else {
      chain_impl();
    }
  }
};

/**
 * Specialization of var multiplication for arithmetic and `var_value`
 */
template <typename VariVal, typename Arith, typename Vari>
class multiply_vari<VariVal, Arith, Vari, require_vt_arithmetic<Arith>> final
    : public op_vari<VariVal, Arith, Vari*> {
  using op_vari<VariVal, Arith, Vari*>::ad;
  using op_vari<VariVal, Arith, Vari*>::bvi;

 public:
  multiply_vari(const Arith& a, Vari* bvi)
      : op_vari<VariVal, Arith, Vari*>(a * bvi->val_, a, bvi) {}

  template <typename T1 = Arith, typename T2 = Vari,
            require_arithmetic_t<T1>* = nullptr,
            require_vari_vt<std::is_arithmetic, T2>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += this->adj_ * ad();
  }

  template <typename T1 = Arith, typename T2 = Vari,
            require_vt_arithmetic<T1>* = nullptr,
            require_vari_vt<is_eigen, T2>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += (this->adj_ * ad()).transpose();
  }
  // NOTE: THIS IS WRONG
  template <typename T1 = Arith, typename T2 = Vari,
            require_eigen_t<T1>* = nullptr,
            require_vari_vt<std::is_arithmetic, T2>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += this->adj_.sum();
  }

  void chain() {
    if (unlikely(is_any_nan(bvi()->val_, ad()))) {
      fill(bvi()->adj_, NOT_A_NUMBER);
    } else {
      chain_impl();
    }
  }
};

}  // namespace internal

template <typename T> 
using require_scalar_t = require_t<std::is_same<std::decay_t<T>, scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_matrix_t = require_t<disjunction<is_eigen<T>,
					       std::is_same<std::decay_t<T>, var_value<Eigen::MatrixXd>>,
					       std::is_same<std::decay_t<T>, var_value<Eigen::VectorXd>>,
					       std::is_same<std::decay_t<T>, var_value<Eigen::RowVectorXd>>>>;

struct OpMultiplyScalarScalar {
  double a_;
  double b_;

  template <std::size_t size>
  double operator()(const std::array<bool, size>& needs_adj,
		    double a,
		    double b) {
    a_ = a;
    b_ = b;

    return a * b;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 double adj) {
    return std::make_tuple(adj * b_, adj * a_);
  }
};

struct OpMultiplyMatrixScalar {
  int N_;
  int M_;
  double* x_mem_;
  double b_;

  template <std::size_t size, typename Derived>
  Eigen::MatrixXd operator()(const std::array<bool, size>& needs_adj,
			     const Eigen::MatrixBase<Derived>& x,
			     double b) {
    N_ = x.rows();
    M_ = x.cols();

    if(needs_adj[1]) {
      x_mem_
        = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(N_ * M_);

      for (int n = 0; n < N_ * M_; ++n) {
	x_mem_[n] = x(n);
      }
    }

    b_ = b;

    return x * b;
  }

  template <std::size_t size, int R, int C>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::MatrixXd& adj) {
    Eigen::MatrixXd adja;
    double adjb = 0.0;

    if(needs_adj[0]) {
      adja.resize(N_, M_);
      adja = adj * b_;
    }
    
    if(needs_adj[1]) {
      Eigen::Map<Eigen::MatrixXd> x(x_mem_, N_, M_);
      adjb = x.dot(adj);
    }

    return std::make_tuple(adja, adjb);
  }
};

struct OpMultiplyMatrixMatrix {
  int N1_;
  int M1_;
  int N2_;
  int M2_;
  double* A_mem_;
  double* B_mem_;

  template <std::size_t size, typename Derived1, typename Derived2>
  Eigen::MatrixXd operator()(const std::array<bool, size>& needs_adj,
			     const Eigen::MatrixBase<Derived1>& A,
			     const Eigen::MatrixBase<Derived2>& B) {
    N1_ = A.rows();
    M1_ = A.cols();
    N2_ = B.rows();
    M2_ = B.cols();

    if(needs_adj[0]) {
      B_mem_
        = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(N2_ * M2_);

      for (int n = 0; n < N2_ * M2_; ++n) {
	B_mem_[n] = B(n);
      }
    }

    if(needs_adj[1]) {
      A_mem_
        = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(N1_ * M1_);

      for (int n = 0; n < N1_ * M1_; ++n) {
	A_mem_[n] = A(n);
      }
    }

    return A * B;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::MatrixXd& adj) {
    Eigen::MatrixXd adjA;
    Eigen::MatrixXd adjB;

    if(needs_adj[0]) {
      Eigen::Map<Eigen::MatrixXd> B(B_mem_, N2_, M2_);
      adjA = adj * B.transpose();
    }
    
    if(needs_adj[1]) {
      Eigen::Map<Eigen::MatrixXd> A(A_mem_, N1_, M1_);
      adjB = A.transpose() * adj;
    }

    return std::make_tuple(adjA, adjB);
  }
};

template <typename T1, typename T2,
	  require_scalar_t<T1>...,
	  require_scalar_t<T2>...,
	  require_any_var_value_t<T1, T2>...>
inline auto operator*(const T1& a, const T2& b) {
  return adj_jac_apply<OpMultiplyScalarScalar>(a, b);
}

template <typename T1, typename T2,
	  require_matrix_t<T1>...,
	  require_scalar_t<T2>...,
	  require_any_var_value_t<T1, T2>...>
inline auto operator*(const T1& a, const T2& b) {
  return adj_jac_apply<OpMultiplyMatrixScalar>(a, b);
}

template <typename T1, typename T2,
	  require_scalar_t<T1>...,
	  require_matrix_t<T2>...,
	  require_any_var_value_t<T1, T2>...>
inline auto operator*(const T1& a, const T2& b) {
  return b * a;
}
 
template <typename T1, typename T2,
	  require_matrix_t<T1>...,
	  require_matrix_t<T2>...,
	  require_any_var_value_t<T1, T2>...>
inline auto operator*(const T1& a, const T2& b) {
  return adj_jac_apply<OpMultiplyMatrixMatrix>(a, b);
}


}  // namespace math
}  // namespace stan
#endif
