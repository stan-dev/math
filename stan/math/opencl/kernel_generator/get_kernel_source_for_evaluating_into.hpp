#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/colwise_reduction.hpp>
#include <CL/cl2.hpp>
#include <string>
#include <set>

namespace stan {
namespace math {

template <typename Derived, typename Scalar, typename... Args>
template <typename T_lhs>
std::string
operation_cl<Derived, Scalar, Args...>::get_kernel_source_for_evaluating_into(
    const T_lhs& lhs) const {
  static_assert(
      is_valid_expression<T_lhs>::value,
      "operation_cl::get_kernel_source_for_evaluating_into: left hand "
      "side is not a valid expression!");
  auto lhs_expression = as_operation_cl(lhs);
  std::set<const operation_cl_base*> generated;
  name_generator ng;
  kernel_parts parts = derived().get_whole_kernel_parts(generated, ng, "i", "j",
                                                        lhs_expression);
  std::string src = "kernel void calculate(" + parts.args +
                    "const int rows, const int cols){\n"
                    "const int i = get_global_id(0);\n"
                    "const int j = get_global_id(1);\n"
                    "const int lid_i = get_local_id(0);\n"
                    "const int lsize_i = get_local_size(0);\n"
                    "const int wg_id_i = get_group_id(0);\n"
                    "const int n_groups_i = get_num_groups(0);\n"
                    + parts.initialization +
                    + "if(i < rows && j < cols){\n" + parts.body + "}\n"
                    + parts.reduction + "}";
  return src;
}

}  // namespace math
}  // namespace stan

#endif
#endif
