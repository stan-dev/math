#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
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
  kernel_parts parts = derived().get_kernel_parts(generated, ng, "i", "j");
  kernel_parts out_parts
      = lhs_expression.get_kernel_parts_lhs(generated, ng, "i", "j");
  std::string src = "kernel void calculate(" + parts.args +
                    out_parts.args.substr(0, out_parts.args.size() - 2) +
                    "){\n"
                    "int i = get_global_id(0);\n"
                    "int j = get_global_id(1);\n"
                    + parts.body +
                    out_parts.body + " = " + var_name + ";\n}";
  return src;
}

}  // namespace math
}  // namespace stan

#endif
#endif
