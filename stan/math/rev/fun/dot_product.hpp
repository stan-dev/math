#ifndef STAN_MATH_REV_FUN_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_DOT_PRODUCT_HPP

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



#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

struct OpDotProductEigen {
  int N_;
  double* a_mem_;
  double* b_mem_;

  template <std::size_t size,
	    typename Derived1,
	    typename Derived2>
  double operator()(const std::array<bool, size>& needs_adj,
		    const Eigen::MatrixBase<Derived1>& a,
		    const Eigen::MatrixBase<Derived2>& b) {
    N_ = a.size();

    double ret = 0.0;

    if(needs_adj[0])
      b_mem_
        = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(N_);
    
    if(needs_adj[1])
      a_mem_
        = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(N_);

    for (int n = 0; n < N_; ++n) {
      ret += a(n) * b(n);

      if(needs_adj[0])
	b_mem_[n] = b(n);

      if(needs_adj[1])
	a_mem_[n] = a(n);
    }

    return ret;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const double& adj) {
    Eigen::VectorXd adja;
    Eigen::VectorXd adjb;

    if(needs_adj[0]) {
      Eigen::Map<Eigen::VectorXd> b(b_mem_, N_);
      adja.resize(N_);
      adja = adj * b;
    }
    
    if(needs_adj[1]) {
      Eigen::Map<Eigen::VectorXd> a(a_mem_, N_);
      adjb.resize(N_);
      adjb = adj * a;
    }

    return std::make_tuple(adja, adjb);
  }
};

}  // namespace internal

/**
 * Returns the dot product.
 *
 * @tparam T1 type of the first vector
 * @tparam T2 type of the second vector
 *
 * @param[in] v1 First row or column vector.
 * @param[in] v2 Second row or column vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error if length of v1 is not equal to length of v2.
 */
template <typename Vec1, typename Vec2,
	  require_t<
	    disjunction<conjunction<
			  require_vector2<Vec1>,
			  require_vector2<Vec2>>,
			conjunction<
			  require_row_vector2<Vec1>,
			  require_row_vector2<Vec2>>>>* = nullptr,
          require_any_var2_t<Vec1, Vec2>* = nullptr>
inline var dot_product(const Vec1& v1, const Vec2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  return adj_jac_apply<internal::OpDotProductEigen>(v1, v2);
}

/**
 * Returns the dot product.
 *
 * @tparam T1 type of elements in the first vector
 * @tparam T2 type of elements in the second vector
 *
 * @param[in] v1 First array.
 * @param[in] v2 Second array.
 * @param[in] length Length of both arrays.
 * @return Dot product of the arrays.
 */
template <typename T1, typename T2,
	  require_any_var_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> dot_product(const T1* v1, const T2* v2,
                                         size_t length) {
  return adj_jac_apply<internal::OpDotProductEigen>
    (Eigen::Map<Eigen::VectorXd>(v1, length),
     Eigen::Map<Eigen::VectorXd>(v2, length));
}
  
/**
 * Returns the dot product.
 *
 * @tparam T1 type of elements in the first vector
 * @tparam T2 type of elements in the second vector
 *
 * @param[in] v1 First vector.
 * @param[in] v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error if sizes of v1 and v2 do not match.
 */
template <typename T1,
	  typename T2,
	  require_any_var_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> dot_product(const std::vector<T1>& v1,
                                         const std::vector<T2>& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  return adj_jac_apply<internal::OpDotProductEigen>
    (Eigen::Map<Eigen::VectorXd>(v1.data(), v1.size()),
     Eigen::Map<Eigen::VectorXd>(v2.data(), v2.size()));
}

}  // namespace math
}  // namespace stan
#endif
