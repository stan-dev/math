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
  std::string src;
  if (Derived::require_specific_local_size) {
    src =
        "kernel void calculate(" + parts.args +
        "const int rows, const int cols){\n"
        "const int gid_i = get_global_id(0);\n"
        "const int lid_i = get_local_id(0);\n"
        "const int lsize_i = get_local_size(0);\n"
        "const int wg_id_i = get_group_id(0);\n"
        "const int wg_id_j = get_group_id(1);\n"
        "const int n_groups_i = get_num_groups(0);\n"
        "const int blocks_rows = (rows + lsize_i - 1) / lsize_i;\n"
        "const int blocks_cols = (cols + lsize_i - 1) / lsize_i;\n"
        "const int i0 = lsize_i * wg_id_i;\n"
        "const int i = i0 + lid_i;\n"
        "const int j0 = lsize_i * wg_id_j;\n"
        "for(int lid_j = 0; lid_j < min(cols - j0, lsize_i); lid_j++){\n"
        "const int j = j0 + lid_j;\n"
        + parts.initialization +
        "if(i < rows){\n"
        + parts.body +
        "}\n"
        + parts.reduction +
        "}\n"
        "}\n";
  } else {
    src =
        "kernel void calculate(" + parts.args.substr(0, parts.args.size() - 2) +
        "){\n"
        "int i = get_global_id(0);\n"
        "int j = get_global_id(1);\n"
        + parts.initialization
        + parts.body
        + parts.reduction +
        "}\n";
  }
  return src;
}

}  // namespace math
}  // namespace stan

#endif
#endif
