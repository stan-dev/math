#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <CL/cl2.hpp>
#include <string>
#include <set>

namespace stan {
namespace math {

template <typename Derived, typename Scalar, typename... Args>
template <typename T_lhs>
void operation_cl<Derived, Scalar, Args...>::evaluate_into(
    const T_lhs& lhs) const {
  static_assert(
      is_valid_expression<T_lhs>::value,
      "operation_cl::evaluate_into: left hand side is not a valid expression!");
  using cache = operation_cl<Derived, Scalar, Args...>::cache<T_lhs>;
  const auto& lhs_expression = as_operation_cl(lhs);

  int n_rows = derived().rows();
  int n_cols = derived().cols();
  const char* function = "evaluate_into";
  check_positive(function, "number of rows", n_rows);
  check_positive(function, "number of columns", n_cols);
  check_size_match(function, "Rows of ", "*this", n_rows, "rows of ",
                   "lhs_expression", lhs_expression.rows());
  check_size_match(function, "Columns of ", "*this", n_cols, "columns of ",
                   "lhs_expression", lhs_expression.cols());
  try {
    if (cache::kernel() == nullptr) {
      cache::source = get_kernel_source_for_evaluating_into(lhs);
      auto opts = opencl_context.base_opts();
      cache::kernel = opencl_kernels::compile_kernel(
          "calculate", {view_kernel_helpers, cache::source.c_str()}, opts);
    }
    int arg_num = 0;
    std::set<const operation_cl_base*> generated;
    derived().set_args(generated, cache::kernel, arg_num);
    lhs_expression.set_args(generated, cache::kernel, arg_num);

    cl::Event e;
    opencl_context.queue().enqueueNDRangeKernel(cache::kernel, cl::NullRange,
                                                cl::NDRange(n_rows, n_cols),
                                                cl::NullRange, nullptr, &e);
    derived().add_read_event(e);
    lhs_expression.add_write_event(e);
  } catch (cl::Error e) {
    check_opencl_error("operation_cl.evaluate_into", e);
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
