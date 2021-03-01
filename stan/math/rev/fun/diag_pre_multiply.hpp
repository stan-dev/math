#ifndef STAN_MATH_REV_FUN_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_DIAG_PRE_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/multiply.hpp>

#include <iostream>


namespace stan {
namespace math {

/**
 * Return the elementwise multiplication of the specified
 * matrices.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression
 * @param m2 Second matrix or expression
 * @return Elementwise product of matrices.
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2>* = nullptr>
/* template <typename Mat1, typename Mat2, 
					require_eigen_vector_t<Mat1>* = nullptr,
          require_eigen_t<Mat2>* = nullptr> */
// template <typename Mat1, typename Mat2>
auto diag_pre_multiply(const Mat1& m1, const Mat2& m2) {
	std::cout << "I am using rev." << std::endl;
  // check_matching_dims("elt_multiply", "m1", m1, "m2", m2);
	
	
	
	check_size_match("diag_pre_multiply", "m1.size()", m1.size(), "m2.rows()",
                   m2.rows());
  using inner_ret_type = decltype(value_of(m1).asDiagonal() * value_of(m2));
  using ret_type = return_var_matrix_t<inner_ret_type, Mat1, Mat2>;
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.val().asDiagonal() * arena_m2.val());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      for (int i = 0; i < arena_m2.cols(); ++i) {
        for (int j = 0; j < arena_m2.rows(); ++j) {
          const auto ret_adj = ret.adj().coeffRef(i, j);
          arena_m1.adj().coeffRef(i) += arena_m2.val().coeff(i, j) * ret_adj;
          arena_m2.adj().coeffRef(i, j) += arena_m1.val().coeff(i) * ret_adj;
        }
      }
    });
    return ret_type(ret);
  } else if (!is_constant<Mat1>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, Mat2>> arena_m2 = value_of(m2);
    arena_t<ret_type> ret(arena_m1.val().asDiagonal() * arena_m2);
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
			for (int i = 0; i < arena_m1.size(); ++i) {
				for (int j = 0; j < arena_m1.size(); ++j) {
					arena_m1.adj().coeffRef(i) += arena_m2.val().coeff(i, j) *
						ret.adj().coeffRef(i, j);
				}
			}
    });
    return ret_type(ret);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
		arena_t<ret_type> ret(arena_m1.asDiagonal() * arena_m2.val());
		std::cout << "Initialized adj():" << std::endl << arena_m2.adj()
			<< std::endl << "-------------------" << std::endl;
		std::cout << "Initialized adj().coeffRef:" << std::endl << arena_m2.adj().coeffRef(0, 0)
			<< std::endl << "-------------------" << std::endl;
		std::cout << "arena_m1.val():" << std::endl << arena_m1.val()
			<< std::endl << "-------------------" << std::endl;
		std::cout << "ret.adj():" << std::endl << ret.adj()
			<< std::endl << "-------------------" << std::endl;
		std::cout << "m1.size():" << std::endl << arena_m1.size()
			<< std::endl << "-------------------" << std::endl;
		std::cout << "ret:" << std::endl << ret
			<< std::endl << "-------------------" << std::endl;
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
			std::cout << "hh" << std::endl;
			/* for (int i = 0; i < arena_m1.size(); ++i) {
				for (int j = 0; j < arena_m1.size(); ++j) {
					std::cout << "oo" << std::endl;
					// arena_m2.adj().coeffRef(i,j) += 2;
					//arena_m2.adj().coeffRef(i,j) += arena_m1.val().coeff(i); *
					//	ret.adj().coeffRef(i, j);
				}
			} */
    });
		std::cout << "Works here." << std::endl;
    return ret_type(ret);
  }
}

}  // namespace math
}  // namespace stan


#endif