#ifndef STAN_MATH_PRIM_META_FWD_DECL_HPP
#define STAN_MATH_PRIM_META_FWD_DECL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/return_type.hpp>

namespace stan {

class operation_cl_base;
class operation_cl_lhs_base;

/**
 * scalar_seq_view provides a uniform sequence-like wrapper around either a
 * scalar or a sequence of scalars.
 *
 * @tparam C the container type; will be the scalar type if wrapping a scalar
 * @tparam T the scalar type
 */
template <typename C, typename = void>
class scalar_seq_view;

/**
 * This class provides a low-cost wrapper for situations where you either need
 * an Eigen Vector or RowVector or a std::vector of them and you want to be
 * agnostic between those two options. This is similar to scalar_seq_view but
 * instead of being a sequence-like view over a scalar or seq of scalars, it's
 * a sequence-like view over a Vector or seq of Vectors. Notably this version
 * only allows std::vectors as the outer container type, since we would have
 * difficulty figuring out which contained type was the container otherwise.
 *
 * @tparam T the wrapped type, either a Vector or std::vector of them.
 */
template <typename T, typename = void>
class vector_seq_view;

namespace math {
namespace internal {

template <typename T, typename S, typename Enable>
class empty_broadcast_array;

/** \ingroup type_trait
 * \callergraph
 * An edge holds both the operands and its associated
 * partial derivatives. They're held together in the
 * same class because then we can keep the templating logic that
 * specializes on type of operand in one place.
 *
 * This is the base template class that ends up getting instantiated
 * for arithmetic primitives (doubles and ints).
 *
 * NB: since ops_partials_edge.partials_ and ops_partials_edge.partials_vec
 * are sometimes represented internally as a broadcast_array, we need to take
 * care with assignments to them. Indeed, we can assign any right hand side
 * which allows for indexing to a broadcast_array. The resulting behaviour is
 * that the entry for the first index is what gets assigned. The most common
 * use-case should be where the rhs is some container of length 1.
 *
 * @tparam ViewElt the type we expect to be at partials_[i]
 * @tparam Op the type of the operand
 */
template <typename ViewElt, typename Op, typename Enable = void>
class ops_partials_edge;

template <typename ReturnType, typename Enable, typename... Ops>
class partials_propagator;

}

/** \ingroup type_trait
 * Primary template class for the metaprogram to compute the index
 * type of a container.
 *
 * Only the specializations have behavior that can be used, and
 * all implement a typedef <code>type</code> for the type of the
 * index given container <code>T</code>.
 *
 * @tparam T type of container.
 */
template <typename T, typename = void>
struct index_type {
  static_assert(1, "index_type not specialized for this type");
  using type = void;
};

template <typename T>
using index_type_t = typename index_type<T>::type;


/**
 * Base template class for vectorization of unary scalar functions
 * defined by a template class <code>F</code> to a scalar,
 * standard library vector, or Eigen dense matrix expression
 * template.
 *
 * <p>Specializations of this base define applications to scalars
 * (primitive or autodiff in the corresponding autodiff library
 * directories) or to standard library vectors of vectorizable
 * types (primitives, Eigen dense matrix expressions, or further
 * standard vectors).
 *
 * <p>Each specialization must define the typedef
 * <code>return_t</code> for the vectorized return type and the
 * function <code>apply</code> which defines the vectorization or
 * base application of the function defined statically by the
 * class F.  The function definition class F defines a static
 * function <code>fun()</code>, which defines the function's
 * behavior on scalars.
 *
 * @tparam F Type of function to apply.
 * @tparam T Type of argument to which function is applied.
 */
template <typename F, typename T, typename Enable = void>
struct apply_scalar_unary;

template <bool ArithmeticArguments, typename F, typename T_y0, typename... Args>
struct coupled_ode_system_impl;

// Forward declaration to allow specializations
template <typename T, typename Enable = void>
struct apply_vector_unary;


/**
 * LDLT_factor is a structure that holds a matrix of type T and the
 * LDLT of its values.
 *
 * @tparam T type of elements in the matrix
 */
template <typename T, typename Enable = void>
class LDLT_factor;


template <typename F, typename T_shared_param, typename T_job_param>
class map_rect_reduce;

template <typename Op1 = double, typename Op2 = double, typename Op3 = double,
          typename Op4 = double, typename Op5 = double, typename Op6 = double,
          typename Op7 = double, typename Op8 = double,
          typename T_return_type
          = return_type_t<Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8>>
class operands_and_partials;  // Forward declaration

/**
 * \defgroup eigen_expressions Eigen expressions
 */
/**
 * \ingroup eigen_expressions
 * \defgroup returning_expressions Returning expressions
 */
/**
 * \ingroup returning_expressions
 *
 * Operations in Eigen Expressions hold their arguments either by value or by
 * reference. Which one is chosen depends on type of the argument. Other
 * operations are held by value. "Heavy" objects that can hold data themselves,
 * such as `Eigen::Matrix` or `Eigen::Ref` are instead held by reference. THis
 * is the only criterion - holding rvalue arguments by value is not supported,
 * so we can not use perfect forwarding.
 *
 * When returning an expression from function we have to be careful that any
 * arguments in this expression that are held by reference do not go out of
 * scope. For instance, consider the function:
 *
 * ```
 * template<typename T>
 * auto f(const T& x){
 *   const Eigen::Ref<const Eigen::VectorXd>& x_ref = x;
 *   return x_ref.cwiseProduct(x_ref);
 * }
 * ```
 * And the code calling it:
 * ```
 * Eigen::MatrixXd test_mat(2,2);
 * test_mat << 5, 5, 5, 5;
 * VectorXd X  = f(test_mat.diagonal());
 * ```
 * This function will return back a `CwiseBinaryOp` Eigen expression, which is
 * then evaluated out of the function scope when assigned to `X`. The expression
 * references `x_ref`, which was created withing function and destroyed when
 * the function returned. The returned expression is evaluated after the
 * function returned, so its evaluation references a matrix that was already
 * deleted. In other words the returned expression contains a dangling
 * reference.
 *
 * So a function returning an expression referencing local matrices or
 * matrices that were rvalue reference arguments to the function will not work.
 *
 * A workarount to this issue is allocating and constructing or moving such
 * objects to heap. `Holder` object is a no-op operation that can also take
 * pointers to such objects and release them when it goes out of scope. It can
 * be created either by directly supplying pointers to such objects to `holder`
 * function or by forwarding function arguments and moving local variables to
 * `make_holder`, which will move any rvalues to heap first.
 */
template <class ArgType, typename... Ptrs>
class Holder;

template <typename T>
class arena_matrix_cl;

/**
 * Non-templated base class for `matrix_cl` simplifies checking if something is
 * matrix_cl.
 */
class matrix_cl_base;

/**
 * Type for sizes and indexes in an Eigen matrix with double elements.
 */
using size_type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Index;

/**
 * Type for matrix of double values.
 */
using matrix_d = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * Type for (column) vector of double values.
 */
using vector_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;

/**
 * Type for (row) vector of double values.
 */
using row_vector_d = Eigen::Matrix<double, 1, Eigen::Dynamic>;


}
}

#endif
