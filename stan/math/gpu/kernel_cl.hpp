#ifndef STAN_MATH_GPU_KERNEL_CL_HPP
#define STAN_MATH_GPU_KERNEL_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/constants.hpp>
#include <stan/math/gpu/kernels/helpers.hpp>
#include <CL/cl.hpp>
#include <string>
#include <algorithm>
#include <map>
#include <vector>

// Used for importing the opencl kernels at compile time.
// There has been much discussion about the best ways to do this:
// https://github.com/bstatcomp/math/pull/7
// and https://github.com/stan-dev/math/pull/966
#ifndef STRINGIFY
#define STRINGIFY(src) #src
#endif

namespace stan {
namespace math {
namespace {

// Holds Default parameter values for each Kernel.
typedef std::map<const char*, int> map_base_opts;
map_base_opts base_opts
    = {{"LOWER", static_cast<int>(TriangularViewGPU::Lower)},
       {"UPPER", static_cast<int>(TriangularViewGPU::Upper)},
       {"ENTIRE", static_cast<int>(TriangularViewGPU::Entire)},
       {"UPPER_TO_LOWER", static_cast<int>(TriangularMapGPU::UpperToLower)},
       {"LOWER_TO_UPPER", static_cast<int>(TriangularMapGPU::LowerToUpper)},
       {"WORK_PER_THREAD_MULT", 8},
       {"WORKGROUP_SIZE_MULT", 32},
       {"WORKGROUP_SIZE_MULT_SELF_TRANS", 32},
       {"WORK_PER_THREAD_MULT_SELF_TRANS", 4}};

auto compile_kernel(const char* name, const char* source) {
  std::string kernel_opts = "";
  for (auto&& comp_opts : base_opts) {
    kernel_opts += std::string(" -D") + comp_opts.first + "="
                   + std::to_string(comp_opts.second);
  }
  std::string kernel_source(opencl_kernels::helpers);
  kernel_source.append(source);
  cl::Program program;
  try {
    cl::Program::Sources src(1, std::make_pair(kernel_source.c_str(),
                                               strlen(kernel_source.c_str())));
    program = cl::Program(opencl_context.context(), src);
    program.build({opencl_context.device()}, kernel_opts.c_str());

    return cl::Kernel(program, name);
  } catch (const cl::Error& e) {
    std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(
        opencl_context.device()[0]);
    std::cerr << "Build log :" << std::endl << buildlog << std::endl;
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

  auto operator()() const { return cl::make_kernel<Args...>(kernel_); }
};

template <typename... Args>
struct global_range_kernel {
  const kernel_functor<Args...> make_functor;
  global_range_kernel(const char* name, const char* source)
      : make_functor(name, source) {}
  auto operator()(cl::NDRange thread_global_size, Args... args) const {
    auto f = make_functor();
    cl::EnqueueArgs eargs(opencl_context.queue(), thread_global_size);
    f(eargs, args...).wait();
  }
};

template <typename... Args>
struct local_range_kernel {
  const kernel_functor<Args...> make_functor;
  local_range_kernel(const char* name, const char* source)
      : make_functor(name, source) {}
  auto operator()(cl::NDRange thread_global_size, cl::NDRange thread_block_size,
                  Args... args) const {
    auto f = make_functor();
    cl::EnqueueArgs eargs(opencl_context.queue(), thread_global_size,
                          thread_block_size);
    f(eargs, args...).wait();
  }
};

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
