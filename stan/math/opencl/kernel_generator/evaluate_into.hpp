#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <stan/math/opencl/kernel_generator/as_operation.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <cl.hpp>
#include <string>
#include <set>

namespace stan {
namespace math {

template <typename Derived, typename ReturnScalar>
template <typename T_lhs>
void operation<Derived, ReturnScalar>::evaluate_into(const T_lhs& lhs) const {
  using enable = enable_if_all_valid_expressions<T_lhs>;
  using cache = operation<Derived, ReturnScalar>::cache<T_lhs>;
  const auto& lhs_expression = as_operation(lhs);

  int n_rows = derived().rows();
  int n_cols = derived().cols();
  const char* function = "evaluate_into";
  if (n_rows != dynamic) {
    check_size_match(function, "Rows of ", "*this", n_rows, "rows of ",
                     "lhs_expression", lhs_expression.rows());
  }
  if (n_cols != dynamic) {
    check_size_match(function, "Columns of ", "*this", n_cols, "columns of ",
                     "lhs_expression", lhs_expression.cols());
  }
  try {
    if (cache::kernel() == NULL) {
      std::string src = get_kernel_source_for_evaluating_into(lhs);
      auto opts = opencl_context.base_opts();
      cache::kernel = opencl_kernels::compile_kernel(
          "calculate", {view_kernel_helpers, src.c_str()}, opts);
    }
    int arg_num = 0;
    std::set<const void*> generated;
    derived().set_args(generated, cache::kernel, arg_num);
    lhs_expression.set_args(generated, cache::kernel, arg_num);

    cl::Event e;
    opencl_context.queue().enqueueNDRangeKernel(cache::kernel, cl::NullRange,
                                                cl::NDRange(n_rows, n_cols),
                                                cl::NullRange, nullptr, &e);
    derived().add_event(e);
    lhs_expression.add_write_event(e);
  } catch (cl::Error e) {
    check_opencl_error("operation.evaluate_into", e);
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
