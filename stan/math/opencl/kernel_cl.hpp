#ifndef STAN_MATH_OPENCL_KERNEL_CL_HPP
#define STAN_MATH_OPENCL_KERNEL_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/kernels/helpers.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/arr/fun/vec_concat.hpp>
#include <CL/cl.hpp>
#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <utility>

// Used for importing the OpenCL kernels at compile time.
// There has been much discussion about the best ways to do this:
// https://github.com/bstatcomp/math/pull/7
// and https://github.com/stan-dev/math/pull/966
#ifndef STRINGIFY
#define STRINGIFY(src) #src
#endif

namespace stan {
namespace math {
namespace opencl_kernels {
namespace internal {
/**
 * Extracts the kernel's arguments, used in the global and local kernel
 * constructor.
 * @tparam For this general template the function will just return back the
 * value passed in.
 * @param t The type that will be returned.
 * @return the input t.
 */
template <typename T>
inline const T& get_kernel_args(const T& t) {
  return t;
}

inline const cl::Buffer& get_kernel_args(const stan::math::matrix_cl& m) {
  return m.buffer();
}

template <typename T>
inline void assign_event(const cl::Event&, to_const_matrix_cl_t<T>&) {}

template <>
inline void assign_event<in_buffer>(const cl::Event& e,
                                    const stan::math::matrix_cl& m) {
  m.add_read_event(e);
}

template <>
inline void assign_event<out_buffer>(const cl::Event& e,
                                     const stan::math::matrix_cl& m) {
  m.add_write_event(e);
}

template <>
inline void assign_event<in_out_buffer>(const cl::Event& e,
                                        const stan::math::matrix_cl& m) {
  m.add_read_write_event(e);
}

template <typename T,
          typename std::enable_if_t<std::is_same<T, cl::Event>::value, int> = 0>
inline void assign_events(const T&) {}

/**
 * Adds the event to any matrices in the arguments in the event vector specified
 * by the buffer directionality.
 * @tparam Arg Arguments given during kernel creation that specify the kernel
 * signature.
 * @tparam Args Arguments given during kernel creation that specify the kernel
 * signature.
 * @param new_event The cl::Event generated involving the arguments.
 * @param m Arguments to the kernel that may be matrices or not. Non-matrices
 * ignored.
 * @param args Arguments to the kernel that may be matrices or not. Non-matrices
 * ignored.
 */
template <typename Arg, typename... Args>
inline void assign_events(const cl::Event& new_event,
                          to_const_matrix_cl_t<Arg>& m,
                          to_const_matrix_cl_t<Args>&... args) {
  assign_event<Arg>(new_event, m);
  assign_events<Args...>(new_event, args...);
}

template <typename T>
inline const std::vector<cl::Event> select_events(to_const_matrix_cl_t<T>& t) {
  return std::vector<cl::Event>();
}

template <>
inline const std::vector<cl::Event> select_events<in_buffer>(
    const stan::math::matrix_cl& m) {
  return m.write_events();
}

template <>
inline const std::vector<cl::Event> select_events<out_buffer>(
    const stan::math::matrix_cl& m) {
  return m.read_write_events();
}

template <>
inline const std::vector<cl::Event> select_events<in_out_buffer>(
    const stan::math::matrix_cl& m) {
  return m.read_write_events();
}

}  // namespace internal

/**
 * Compile an OpenCL kernel.
 *
 * @param name The name for the kernel
 * @param sources A std::vector of strings containing the code for the kernel.
 * @param options The values of macros to be passed at compile time.
 */
inline auto compile_kernel(const char* name,
                           const std::vector<const char*>& sources,
                           std::map<const char*, int>& options) {
  std::string kernel_opts = "";
  for (auto&& comp_opts : options) {
    kernel_opts += std::string(" -D") + comp_opts.first + "="
                   + std::to_string(comp_opts.second);
  }
  std::string kernel_source;
  for (const char* source : sources) {
    kernel_source.append(source);
  }
  cl::Program program;
  try {
    cl::Program::Sources src(1, std::make_pair(kernel_source.c_str(),
                                               strlen(kernel_source.c_str())));
    program = cl::Program(opencl_context.context(), src);
    program.build({opencl_context.device()}, kernel_opts.c_str());

    return cl::Kernel(program, name);
  } catch (const cl::Error& e) {
    // in case of CL_BUILD_PROGRAM_FAILURE, print the build error
    if (e.err() == -11) {
      std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(
          opencl_context.device()[0]);
      system_error("compile_kernel", name, e.err(), buildlog.c_str());
    } else {
      check_opencl_error(name, e);
    }
  }
  return cl::Kernel();  // never reached because check_opencl_error throws
}

/**
 * Functor used for compiling kernels.
 *
 * @tparam Args Parameter pack of all kernel argument types.
 */
template <typename... Args>
class kernel_functor {
 private:
  cl::Kernel kernel_;
  std::map<const char*, int> opts_;

 public:
  /**
   * functor to access the kernel compiler.
   * @param name The name for the kernel.
   * @param sources A std::vector of strings containing the code for the kernel.
   * @param options The values of macros to be passed at compile time.
   */
  kernel_functor(const char* name, const std::vector<const char*>& sources,
                 const std::map<const char*, int>& options) {
    auto base_opts = opencl_context.base_opts();
    for (auto& it : options) {
      if (base_opts[it.first] > it.second) {
        base_opts[it.first] = it.second;
      }
    }
    kernel_ = compile_kernel(name, sources, base_opts);
    opts_ = base_opts;
  }

  auto operator()() const { return cl::make_kernel<Args...>(kernel_); }

  /**
   * @return The options that the kernel was compiled with.
   */
  inline const std::map<const char*, int>& get_opts() const { return opts_; }
};

/**
 * Creates functor for kernels
 *
 * @tparam Args Parameter pack of all kernel argument types.
 */
template <typename... Args>
struct kernel_cl {
  const kernel_functor<internal::to_const_buffer_t<Args>&...> make_functor;
  /**
   * Creates functor for kernels that only need access to defining
   *  the global work size.
   * @param name The name for the kernel
   * @param source A string literal containing the code for the kernel.
   * @param options The values of macros to be passed at compile time.
   */
  kernel_cl(const char* name, const char* source,
            const std::map<const char*, int>& options = {})
      : make_functor(name, {source}, options) {}
  /**
   * Creates functor for kernels that only need access to defining
   *  the global work size.
   * @param name The name for the kernel
   * @param sources A std::vector of strings containing the code for the kernel.
   * @param options The values of macros to be passed at compile time.
   */
  kernel_cl(const char* name, const std::vector<const char*>& sources,
            const std::map<const char*, int>& options = {})
      : make_functor(name, sources, options) {}
  /**
   * Executes a kernel
   * @param global_thread_size The global work size.
   * @param args The arguments to pass to the kernel.
   * @tparam Args Parameter pack of all kernel argument types.
   */
  auto operator()(cl::NDRange global_thread_size,
                  internal::to_const_matrix_cl_t<Args>&... args) const {
    auto f = make_functor();
    const std::vector<cl::Event> kernel_events
        = vec_concat(internal::select_events<Args>(args)...);
    cl::EnqueueArgs eargs(opencl_context.queue(), kernel_events,
                          global_thread_size);
    cl::Event kern_event = f(eargs, internal::get_kernel_args(args)...);
    internal::assign_events<Args...>(kern_event, args...);
    return kern_event;
  }

  /**
   * Executes a kernel
   * @param global_thread_size The global work size.
   * @param thread_block_size The thread block size.
   * @param args The arguments to pass to the kernel.
   * @tparam Args Parameter pack of all kernel argument types.
   */
  auto operator()(cl::NDRange global_thread_size, cl::NDRange thread_block_size,
                  internal::to_const_matrix_cl_t<Args>&... args) const {
    auto f = make_functor();
    const std::vector<cl::Event> kernel_events
        = vec_concat(internal::select_events<Args>(args)...);
    cl::EnqueueArgs eargs(opencl_context.queue(), kernel_events,
                          global_thread_size, thread_block_size);
    cl::Event kern_event = f(eargs, internal::get_kernel_args(args)...);
    internal::assign_events<Args...>(kern_event, args...);
    return kern_event;
  }
};

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
