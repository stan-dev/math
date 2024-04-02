#ifndef STAN_MATH_PRIM_META_HPP
#define STAN_MATH_PRIM_META_HPP

/**
 * \defgroup type_trait Type Traits
 * The type traits in Stan math are a mix of custom traits for detecting
 * value and container types of Eigen matrices, standard vectors, standard
 * complex numbers, and backports of C++17 traits.
 */

/**
 * \defgroup require_meta Available requires<> for overloading.
 *
 * See @ref require_meta_doc for an overview of their use and for how to
 *  build your own custom requires for functions.
 *
 * The modules below group the different possible require constraints
 *  available for overloading functions.
 */

/**
 * \ingroup require_meta
 * \defgroup require_stan_scalar Scalar types
 * `require` type traits for types that are either `arithmetic`, `var`, `fvar`,
 * or `Complex`.
 */

/**
 * \ingroup require_stan_scalar
 * \defgroup require_stan_scalar_real Real types
 * `require` type traits for types that are either `arithmetic`, `var`, or
 * `fvar`.
 */

/**
 * \ingroup require_stan_scalar
 * \defgroup require_stan_scalar_complex Complex types
 * `require` type traits for types that are `Complex<T>`.
 */

/**
 * \ingroup require_meta
 * \defgroup require_std Standard library types and traits
 * `require` type traits that come from the standard library.
 */

/**
 * \ingroup require_meta
 * \defgroup require_eigens_types Eigen
 * `require` type traits to detect Eigen types.
 */

/**
 * \ingroup require_meta
 * \defgroup general_types General Types
 * `require` type traits for general types.
 */

/**
 * \ingroup require_meta
 * \defgroup macro_helpers Require Macro Generators
 * These macros are used on type traits to define the set of `requires`
 */

/**
 * \ingroup require_meta
 * \defgroup require_opencl_types OpenCL
 * `require` type traits to detect types used with OpenCL.
 */
#include <stan/math/prim/meta/compiler_attributes.hpp>
#include <stan/math/prim/meta/ad_promotable.hpp>
#include <stan/math/prim/meta/append_return_type.hpp>
#include <stan/math/prim/meta/base_type.hpp>
#include <stan/math/prim/meta/contains_std_vector.hpp>
#include <stan/math/prim/meta/error_index.hpp>
#include <stan/math/prim/meta/forward_as.hpp>
#include <stan/math/prim/meta/holder.hpp>
#include <stan/math/prim/meta/include_summand.hpp>
#include <stan/math/prim/meta/index_type.hpp>
#include <stan/math/prim/meta/index_apply.hpp>
#include <stan/math/prim/meta/is_autodiff.hpp>
#include <stan/math/prim/meta/is_arena_matrix.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/is_dense_dynamic.hpp>
#include <stan/math/prim/meta/is_double_or_int.hpp>
#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/is_constant.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_container_or_var_matrix.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_eigen_dense_base.hpp>
#include <stan/math/prim/meta/is_eigen_dense_dynamic.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/is_eigen_matrix_base.hpp>
#include <stan/math/prim/meta/is_eigen_sparse_base.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/prim/meta/is_matrix_cl.hpp>
#include <stan/math/prim/meta/is_matrix.hpp>
#include <stan/math/prim/meta/is_plain_type.hpp>
#include <stan/math/prim/meta/is_string_convertible.hpp>
#include <stan/math/prim/meta/is_tuple.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_var_and_matrix_types.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/is_var_dense_dynamic.hpp>
#include <stan/math/prim/meta/is_var_eigen.hpp>
#include <stan/math/prim/meta/is_rev_matrix.hpp>
#include <stan/math/prim/meta/is_vari.hpp>
#include <stan/math/prim/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_stan_scalar_or_eigen.hpp>
#include <stan/math/prim/meta/modify_eigen_options.hpp>
#include <stan/math/prim/meta/partials_return_type.hpp>
#include <stan/math/prim/meta/partials_type.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/prim/meta/possibly_sum.hpp>
#include <stan/math/prim/meta/promote_args.hpp>
#include <stan/math/prim/meta/promote_scalar_type.hpp>
#include <stan/math/prim/meta/ref_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/seq_view.hpp>
#include <stan/math/prim/meta/static_select.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/void_t.hpp>
#include <stan/math/prim/meta/StdVectorBuilder.hpp>
#include <stan/math/prim/meta/VectorBuilder.hpp>

#endif
