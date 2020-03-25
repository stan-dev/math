#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_HPP
#ifdef STAN_OPENCL

/**
 * \ingroup opencl
 * \defgroup opencl_kernel_generator OpenCL Kernel Generator
 * OpenCL kernel generator is used to combine multiple matrix operations into a
 * single OpenCL kernel. This is much simpler than writing such a kernel by
 * hand.
 *
 * Using one kernel for multiple operations is faster than using one kernel per
 * operation. Each kernel must load its arguments from global GPU memory
 * and store results back. Memory transfers are relatively slow compared to
 * calculations. By combining mutiple operations into a single kernel some
 * memory transfers can be avoided.
 *
 * Kernel generator uses lazy evaluation. Each operation is represented by an
 * operation object. Such an object holds arguments of the operation. When an
 * operation is assigned to a \c matrix_cl, a left-hand-side operation or \c
 * .eval() is called, an operation is evaluated. Arguments to operations can be
 * other operations, scalars and \c matrix_cl objects.
 *
 * A new kernel generator operation class must be derived from \c operation_cl.
 * Optionally, if the operation should support being assigned to, it can be
 * derived from \c operation_cl_lhs instead. In either case parent template
 * arguments should be set to derived type, type of scalar and types of any
 * expression arguements. Member type \c Scalar should be defined as scalar type
 * of the result of operation. Member function \c generate should accept two
 * strings with expressions for indices in kernel and one string with variable
 * name for each argument to the expression. \c generate should return \c
 * kernel_parts struct that contains code for execution of the operation.
 * Member function \c view should return \c matrix_cl_view of the result. Member
 * function \c deep_copy should make a copy of the expression. Arguments that
 * are operations should be copied by calling their \c deep_copy.
 *
 * Following functions can optionally be defined. Defaults are implemented in \c
 * operation_cl:
 * - <code> void modify_argument_indices(std::string& i, std::string& j)
 * </code>: Modifies what indices are passed to argument's \c
 * generate. By default does nothing.
 * - <code> void set_args(std::set<const operation_cl_base*>& generated,
 * cl::Kernel& kernel, int& arg_num) </code>: Sets kernel arguments. By defult
 * only calls \c set_args on arguments.
 * - <code> int rows() </code>: Returns number of rows of the result. By default
 * returns maximum of the arguments' rows.
 * - <code> int cols() </code>: Returns number of columns of the result. By
 * default returns maximum of the arguments' columns.
 * - <code> int thread_rows() </code>: Number of threads required for this
 * operation in rows direction. By default equals to \c rows().
 * - <code> int thread_cols() </code>: Number of threads required for this
 * operation in cols direction. By default equals to \c cols().
 * - <code> int bottom_diagonal() </code>: Index of bottom nonzero diagonal of
 * the result (0 is the diagonal, positive values are superdiagonals, negative
 * values are subdiagonals). By default returns minimum of arguments \c
 * bottom_diagonal.
 * - <code> int top_diagonal() </code>: Index of top nonzero diagonal of the
 * result (0 is the diagonal, positive values are superdiagonals, negative
 * values are subdiagonals). By default returns maximum of arguments \c
 * top_diagonal.

 * If an operation should support being assigned to it should also define the
 * following. Member function \c generate_lhs with same signature as \c generate
 * that returns code when the operation is assigned to. Following functions can
 * be optionally defined. Defaults are in \c operation_cl_lhs.
 * - <code> void set_view(int bottom_diagonal, int top_diagonal, int
 * bottom_zero_diagonal, int top_zero_diagonal) </code>: Sets view of the
 * underlying matrix depending on which are the extreme sub-/super-diagonals
 * written. By default just calls \c set_view on arguments with same arguments.
 * - <code> void check_assign_dimensions(int rows, int cols) </code>: If the
 * operation size can be modified, it should be set to given size. Otherwise it
 * should check that this operation's size matches given size. By default calls
 * `check_assign_dimensions` on arguments with same arguments.
 *
 * A new operation should also have a user-facing function that accepts
 * arguments to the operation and returns the operation object. Arguments should
 * be passed trough function \c as_operation_cl so that they are wrapped in
 * operations if they are not already operations. If the operation defines
 * \c modify_argument_indices this function should make copies of arguments by
 * calling \c .deep_copy() on them.
 */

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>

#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/kernel_generator/unary_function_cl.hpp>
#include <stan/math/opencl/kernel_generator/block.hpp>
#include <stan/math/opencl/kernel_generator/select.hpp>
#include <stan/math/opencl/kernel_generator/rowwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/colwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/transpose.hpp>

#include <stan/math/opencl/kernel_generator/multi_result_kernel.hpp>
#include <stan/math/opencl/kernel_generator/get_kernel_source_for_evaluating_into.hpp>
#include <stan/math/opencl/kernel_generator/evaluate_into.hpp>

#include <stan/math/opencl/kernel_generator/matrix_cl_conversion.hpp>

#endif
#endif
