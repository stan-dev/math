#ifndef STAN_MATH_GPU_KERNEL_CL_HPP
#define STAN_MATH_GPU_KERNEL_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/constants.hpp>
#include <CL/cl.hpp>
#include <string>
#include <algorithm>
#include <map>
#include <vector>

namespace stan {
namespace math {
namespace {

std::string helpers =                        // Helper macros for the kernels.
#include <stan/math/gpu/kernels/helpers.cl>  // NOLINT
    ;                                        // NOLINT
// Holds Default parameter values for each Kernel.
typedef std::map<const char*, int> map_base_opts;
map_base_opts base_opts
    = {{"LOWER", static_cast<int>(TriangularViewGPU::Lower)},
       {"UPPER", static_cast<int>(TriangularViewGPU::Upper)},
       {"ENTIRE", static_cast<int>(TriangularViewGPU::Entire)},
       {"UPPER_TO_LOWER", static_cast<int>(TriangularMapGPU::UpperToLower)},
       {"LOWER_TO_UPPER", static_cast<int>(TriangularMapGPU::LowerToUpper)}};

auto compile_kernel(const char* name, const char* source) {
  std::string kernel_opts = "";
  for (auto&& comp_opts : base_opts) {
    kernel_opts += std::string(" -D") + comp_opts.first + "="
                   + std::to_string(comp_opts.second);
  }
  std::string kernel_source(helpers);
  kernel_source.append(source);
  try {
    cl::Program::Sources src(1, std::make_pair(kernel_source.c_str(),
                                               strlen(kernel_source.c_str())));
    cl::Program program = cl::Program(opencl_context.context(), src);
    program.build({opencl_context.device()}, kernel_opts.c_str());

    return cl::Kernel(program, name);
  } catch (const cl::Error& e) {
    check_opencl_error(name, e);
  }
  return cl::Kernel();  // never reached because check_opencl_error throws
}
}  // namespace

namespace opencl_kernels {

template <typename... Args>
class kernel_functor {
 private:
  cl::Kernel kernel_;

 public:
  kernel_functor(const char* name, const char* source)
      : kernel_(compile_kernel(name, source)) {}

  auto operator()() const { return cl::make_kernel<Args...>(kernel_) };
};

template <typename... Args>
struct range_2d_kernel {
  const kernel_functor<Args...> make_functor;
  range_2d_kernel(const char* name, const char* source)
      : make_functor(name, source) {}
  auto operator()(cl::NDRange thread_block_size, Args... args) const {
    auto f = make_functor();
    cl::EnqueueArgs eargs(opencl_context.queue(), thread_block_size);
    f(eargs, args...).wait();
  };
};

const range_2d_kernel<cl::Buffer, int, int> identity("identity",
#include <stan/math/gpu/kernels/identity_matrix.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, int, int> copy("copy",
#include <stan/math/gpu/kernels/copy_matrix.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, int, int> transpose("transpose",
#include <stan/math/gpu/kernels/transpose_matrix.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int> add("add",
#include <stan/math/gpu/kernels/add_matrix.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, cl::Buffer, int, int> subtract(
    "subtract",
#include <stan/math/gpu/kernels/subtract_matrix.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, int, int, int, int, int, int, int,
                      int, int, int>
    sub_block("sub_block",
#include <stan/math/gpu/kernels/sub_block.cl>  // NOLINT
    );  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, int, int> check_diagonal_zeros(
    "is_zero_on_diagonal",
#include <stan/math/gpu/kernels/check_diagonal_zeros.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, int, int> check_nan("is_nan",
#include <stan/math/gpu/kernels/check_nan.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, int, int, const double>
    check_symmetric("is_symmetric",
#include <stan/math/gpu/kernels/check_symmetric.cl>  // NOLINT
    );  // NOLINT
const range_2d_kernel<cl::Buffer, cl::Buffer, int, int, TriangularViewGPU>
    copy_triangular("copy_triangular",
#include <stan/math/gpu/kernels/copy_triangular_matrix.cl>  // NOLINT
    );  // NOLINT
const range_2d_kernel<cl::Buffer, int, int, TriangularViewGPU> zeros("zeros",
#include <stan/math/gpu/kernels/zeros_matrix.cl>  // NOLINT
);  // NOLINT
const range_2d_kernel<cl::Buffer, int, int, TriangularMapGPU>
    triangular_transpose("triangular_transpose",
#include <stan/math/gpu/kernels/triangular_transpose.cl>  // NOLINT
    );  // NOLINT

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
