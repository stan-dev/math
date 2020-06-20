#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/prim/meta.hpp>

#include <type_traits>

namespace stan {
namespace math {

namespace internal {

struct ColumnsDotProductOp : public AdjJacOp {
  int N_;
  int M_;
  double* a_mem_;
  double* b_mem_;

  template <std::size_t size,
	    typename Derived1,
	    typename Derived2>
  Eigen::RowVectorXd operator()(const std::array<bool, size>& needs_adj,
				const Eigen::MatrixBase<Derived1>& a,
				const Eigen::MatrixBase<Derived2>& b) {
    N_ = a.rows();
    M_ = a.cols();

    if(needs_adj[0]) {
      b_mem_ = allocate_and_save(b);
    }

    if(needs_adj[1]) {
      a_mem_ = allocate_and_save(a);
    }

    Eigen::RowVectorXd out(M_);
    for(size_t m = 0; m < M_; ++m)
      out(m) = a.col(m).dot(b.col(m));
 
    return out;
  }

  template <std::size_t size, typename Derived>
  decltype(auto) multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::MatrixBase<Derived>& adj) {
    auto a = map_matrix(a_mem_, (needs_adj[1]) ? N_ : 0, M_);
    auto b = map_matrix(b_mem_, (needs_adj[0]) ? N_ : 0, M_);

    return std::make_tuple(b * adj.asDiagonal(),
			   a * adj.asDiagonal());
  }
};

}

template <typename Mat1, typename Mat2,
          require_matrix2_t<Mat1>* = nullptr,
          require_matrix2_t<Mat2>* = nullptr,
          require_any_var2_t<Mat1, Mat2>* = nullptr>
inline auto columns_dot_product(const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);

  return adj_jac_apply<internal::ColumnsDotProductOp>(v1, v2);
}

}  // namespace math
}  // namespace stan
#endif
