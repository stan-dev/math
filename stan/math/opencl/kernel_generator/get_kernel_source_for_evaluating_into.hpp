#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <stan/math/opencl/kernel_generator/as_operation.hpp>
#include <cl.hpp>
#include <string>
#include <set>

namespace stan {
namespace math {

template <typename Derived, typename ReturnScalar>
template <typename T_lhs>
std::string
operation<Derived, ReturnScalar>::get_kernel_source_for_evaluating_into(
    const T_lhs& lhs) const {
  auto lhs_expression = as_operation(lhs);
  std::set<const void*> generated;
  name_generator ng;
  kernel_parts parts = derived().generate(generated, ng, "i", "j");
  kernel_parts out_parts = lhs_expression.generate_lhs(generated, ng, "i", "j");
  std::string src = "kernel void calculate(" + parts.args +
                    out_parts.args.substr(0, out_parts.args.size() - 2) +
                    "){\n"
                    "int i = get_global_id(0);\n"
                    "int j = get_global_id(1);\n"
                    + parts.body +
                    out_parts.body + " = " + var_name + ";}";
  return src;
}

}  // namespace math
}  // namespace stan

#endif
#endif
