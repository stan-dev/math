#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_HPP
#ifdef STAN_OPENCL

/**
 * \ingroup opencl
 * \defgroup opencl_kernel_generator OpenCL Kernel Generator
 *
 * The OpenCL kernel generator is used to combine multiple matrix operations
 * into a single OpenCL kernel. This is much simpler than writing
 * multi-operation kernels by hand.
 *
 * Because global GPU memory loads and stores are relativly slow compared to
 * calculations in a kernel, using one kernel for multiple operations is faster
 * than using one kernel per operation.
 *
 * The kernel generator uses lazy evaluation. Each operation is represented by
 * an object derived from `operation_cl`. Such an object holds arguments of the
 * operations as well as meta information needed to generate calculations on the
 * arguments. Arguments to operations can be other operations, scalars
 * or `matrix_cl` objects. An operation is evaluated when either an operation is
 * assigned to a `matrix_cl` or a left-hand-side operation or `.eval()` is
 * called.
 *
 * ## Defining a new kernel generator operation
 *
 * Element-wise functions can be added using one of the macros in
 * `elt_functions.hpp`.
 *
 * New kernel generator classes must satsify the conditions below:
 *
 * 1. The class must be derived from a class inheriting from `operation_cl`.
 * Optionally, if the operation should support being assigned to, it can be
 * derived from a class inheriting `operation_cl_lhs` instead.
 * 2. It's parent template arguments should be set to derived type, type of
 *  scalar and types of any expression arguements.
 * 3. Member type `Scalar` should be defined as scalar type of the result of
 * the operation.
 * 4. Member function `deep_copy` should make a copy of the expression.
 * Arguments that are operations should be copied by calling their `deep_copy`.
 *
 * The following functions can optionally be defined. Defaults are implemented
 * in `operation_cl`:
 * - `kernel_parts generate(const std::string& i, const std::string& j,
 *                            const std::string& var_name_arg)`:
 *     - Generates kernel code for this expression
 *     - Default: Generates code that copies the result of the first argument.
 * - `void modify_argument_indices(std::string& i, std::string& j)`:
 *     - Modifies what indices are passed to argument's `generate()`.
 *     - Default: No-op
 * - `void set_args(std::set<const operation_cl_base*>& generated,
 *        cl::Kernel& kernel, int& arg_num)`:
 *     - Sets additional kernel arguments.
 *     - Default: Calls `set_args()` on arguments.
 * - `int rows()`:
 *     - Returns Number of rows of the result.
 *     - Default: Returns maximum of the arguments' rows.
 * - `int cols()`:
 *     - Returns number of columns of the result.
 *     - Default: Returns maximum of the arguments' columns.
 * - `int thread_rows()`:
 *     - Returns number of threads required for this operation in rows
 * direction.
 *     - Default: returns `rows()`.
 * - `int thread_cols()`:
 *     - Returns number of threads required for this operation in cols
 * direction.
 *     - Default: `cols()`.
 * - `std::pair<int,int> extreme_diagonals()`:
 *     - Returns indices of the extreme bottom and top nonzero diagonals of the
 * result of the operation. Diagonal has index 0, first superdiagonal 1, first
 * subdiagonal -1 and so on. For instance `load_` of a matrix with lower
 * triangular view would return bottom diagonal of `1-rows()` and top diagonal
 * `0`. Returning a more extreme value is also allowed and especially useful for
 * operations that are broadcast so their exact size is not known.
 *     - Default: Bottom diagonal equals to min of bottom diagonals of
 * arguments. Top diagonal equals to max of top diagonals of arguments.
 *
 * The below functions can be optionally defined for operations that support
 * being assigned to. Defaults are in `operation_cl_lhs`.
 * - `kernel_parts generate_lhs(const std::string& i, const std::string& j,
 *                            const std::string& var_name_arg)`:
 *     - Generates kernel code for this expression
 *     - Default: Generates code that copies the result of the first argument.
 * - `void set_view(int bottom_diagonal, int top_diagonal, int
 * bottom_zero_diagonal, int top_zero_diagonal)`:
 *    - Sets view of the underlying `matrix_cl` depending on where the extreme
 * sub-/super-diagonals  are written.
 *    - Default: Calls `set_view` on expression arguments with the same
 * arguments.
 * - `void check_assign_dimensions(int rows, int cols)`:
 *    - If the operation size can be modified, it should be set to the given
 * size. Otherwise it should check that this operation's size matches given
 * size.
 *    - Default: By default calls `check_assign_dimensions` on expression
 * arguments with the same arguments.
 *
 * A new operation should also have a user-facing function that accepts
 * arguments to the operation and returns the operation object. Arguments should
 * be passed trough function `as_operation_cl` so that they are wrapped in
 * operations if they are not operations themselves. If the operation defines
 * `modify_argument_indices` this function should make copies of arguments by
 * calling `deep_copy()` on them internally.
 */

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>

#include <stan/math/opencl/kernel_generator/as_column_vector_or_scalar.hpp>
#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/kernel_generator/constant.hpp>
#include <stan/math/opencl/kernel_generator/append.hpp>
#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/kernel_generator/elt_function_cl.hpp>
#include <stan/math/opencl/kernel_generator/unary_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/block_zero_based.hpp>
#include <stan/math/opencl/kernel_generator/select.hpp>
#include <stan/math/opencl/kernel_generator/rowwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/colwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/reduction_2d.hpp>
#include <stan/math/opencl/kernel_generator/transpose.hpp>
#include <stan/math/opencl/kernel_generator/broadcast.hpp>
#include <stan/math/opencl/kernel_generator/optional_broadcast.hpp>
#include <stan/math/opencl/kernel_generator/diagonal.hpp>
#include <stan/math/opencl/kernel_generator/holder_cl.hpp>
#include <stan/math/opencl/kernel_generator/check_cl.hpp>
#include <stan/math/opencl/kernel_generator/index.hpp>
#include <stan/math/opencl/kernel_generator/indexing.hpp>
#include <stan/math/opencl/kernel_generator/opencl_code.hpp>
#include <stan/math/opencl/kernel_generator/cast.hpp>

#include <stan/math/opencl/kernel_generator/multi_result_kernel.hpp>
#include <stan/math/opencl/kernel_generator/get_kernel_source_for_evaluating_into.hpp>
#include <stan/math/opencl/kernel_generator/evaluate_into.hpp>

#include <stan/math/opencl/kernel_generator/matrix_cl_conversion.hpp>
#include <stan/math/opencl/kernel_generator/compound_assignments.hpp>

#include <stan/math/opencl/kernel_generator/matrix_vector_multiply.hpp>

#endif
#endif
