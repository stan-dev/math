#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP

#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <stan/math/opencl/kernel_generator/as_operation.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <utility>

namespace stan{
namespace math{

template<typename Derived, typename ReturnScalar>
template<typename T_lhs>
void operation<Derived, ReturnScalar>::evaluate_into(T_lhs&& lhs) const {
  using enable = enable_if_all_valid_expressions<T_lhs>;
  auto lhs_expression = as_operation(std::forward<T_lhs>(lhs));

  int n_rows = derived().rows();
  int n_cols = derived().cols();
  const char* function = "evaluate_into";
  if(n_rows!=dynamic) {
    check_size_match(function, "Rows of ", "*this", n_rows, "rows of ", "lhs_expression", lhs_expression.rows());
  }
  if(n_cols!=dynamic) {
    check_size_match(function, "Columns of ", "*this", n_cols, "columns of ", "lhs_expression", lhs_expression.cols());
  }
  name_generator ng;
  std::set<int> generated;
  kernel_parts parts = derived().generate(ng, generated,"i","j");
  kernel_parts out_parts = lhs_expression.generate_lhs(ng, generated,"i","j");
  std::string src = "kernel void calculate(" + parts.args + out_parts.args.substr(0,out_parts.args.size()-2) + "){\n"
                   "int i = get_global_id(0);"
                   "int j = get_global_id(1);\n"
                    + parts.body +
                    out_parts.body + " = " + var_name + ";}";
  try {
    if(kernel_cache.count(src)==0){
      auto opts = opencl_context.base_opts();
      kernel_cache[src] = opencl_kernels::compile_kernel("calculate", {view_kernel_helpers, src.c_str()}, opts);
    }
    cl::Kernel& kernel = kernel_cache[src];
    int arg_num = 0;
    generated.clear();
    derived().set_args(generated,kernel,arg_num);
    lhs_expression.set_args(generated, kernel, arg_num);

    cl::Event e;
    opencl_context.queue().enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(n_rows,n_cols), cl::NullRange, nullptr, &e);
    derived().add_event(e);
    lhs_expression.add_write_event(e);
  }
  catch(cl::Error e){
    check_opencl_error("block.operator=", e);
  }
}

}
}


#endif
