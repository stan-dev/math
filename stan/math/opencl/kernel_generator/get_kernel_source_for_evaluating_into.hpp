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
  std::string src =
      "kernel void calculate(" + parts.args + "const int rows, const int cols){\n"
      "const int gid = get_global_id(0);\n"
      "const int lid = get_local_id(0);\n"
      "const int lsize = get_local_size(0);\n"
      "const int wg_id = get_group_id(0);\n"
      "const int n_groups = get_num_groups(0);\n"

      "const int blocks_rows = (rows + lsize - 1) / lsize;\n"
      "const int blocks_cols = (cols + lsize - 1) / lsize;\n"

      "for (int idx = wg_id; idx < blocks_rows * blocks_cols; idx += n_groups){\n"
      "const int i0 = lsize * (idx % blocks_rows);\n"
      "const int j0 = lsize * (idx / blocks_rows);\n"
      "const int local_rows = min(lsize, rows - i0);\n"
      "const int local_cols = min(lsize, cols - j0);\n"
      "int local_work = (local_rows * local_cols + lsize - 1) / lsize * lsize;\n"

      "for(int idx_local = lid; idx_local < local_work; idx_local += lsize){\n"
      "const int i_local = idx_local % local_rows;\n"
      "const int j_local = idx_local / local_rows;\n"
      "const int idx_local_min = idx_local - lid;\n"
      "const int i_local_min = idx_local_min % local_rows;\n"
      "const int j_local_min = idx_local_min / local_rows;\n"
      "const int idx_local_max = idx_local - lid + lsize - 1;\n"
      "const int i_local_max = idx_local_max % local_rows;\n"
      "const int j_local_max = idx_local_max / local_rows;\n"
      "const int i = i0 + i_local;\n"
      "const int j = j0 + j_local;\n"
      + parts.initialization +
      "if(i < rows && j < cols){\n"
      + parts.body +
      "}\n"
      + parts.reduction +
      "}\n"
      "}\n"
      "}";
  return src;
}

}  // namespace math
}  // namespace stan

#endif
#endif
