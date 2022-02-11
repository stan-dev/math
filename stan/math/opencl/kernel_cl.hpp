#ifndef STAN_MATH_OPENCL_KERNEL_CL_HPP
#define STAN_MATH_OPENCL_KERNEL_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/stringify.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/kernels/helpers.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <CL/opencl.hpp>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <utility>

namespace stan {
namespace math {
namespace opencl_kernels {
namespace internal {

/** \ingroup kernel_executor_opencl
 * Extracts the kernel's arguments, used in the global and local kernel
 * constructor.
 * @tparam For this general template the function will return back the
 * value passed in.
 * @param t The type that will be returned.
 * @return the input t.
 */
template <typename T, require_not_matrix_cl_t<T>* = nullptr>
inline const T& get_kernel_args(const T& t) {
  return t;
}

/** \ingroup kernel_executor_opencl
 * Extracts the kernel's arguments, used in the global and local kernel
 * constructor.
 * @tparam K The type of the \c matrix_cl.
 * @param m The \c matrix with an OpenCL Buffer.
 * @return The OpenCL Buffer.
 */
template <typename K, require_matrix_cl_t<K>* = nullptr>
inline const cl::Buffer& get_kernel_args(const K& m) {
  return m.buffer();
}

/** \ingroup kernel_executor_opencl
 * Assigns the event to a \c matrix_cl.
 * @tparam T The type to be assigned, if not a matrix_cl this function
 * will do nothing.
 * @tparam K The type of the \c matrix_cl.
 * @param e The event to be assigned.
 */
template <typename T, require_not_matrix_cl_t<T>* = nullptr>
inline void assign_event(const cl::Event& e, const T&) {}

template <typename T, typename K, require_matrix_cl_t<K>* = nullptr,
          require_same_t<T, in_buffer>* = nullptr>
inline void assign_event(const cl::Event& e, const K& m) {
  m.add_read_event(e);
}
template <typename T, typename K, require_matrix_cl_t<K>* = nullptr,
          require_same_t<T, out_buffer>* = nullptr>
inline void assign_event(const cl::Event& e, K& m) {
  m.add_write_event(e);
}
template <typename T, typename K, require_matrix_cl_t<K>* = nullptr,
          require_same_t<T, in_out_buffer>* = nullptr>
inline void assign_event(const cl::Event& e, K& m) {
  m.add_read_write_event(e);
}

/** \ingroup kernel_executor_opencl
 * Adds the event to any \c matrix_cls in the arguments depending on whether
 * they are \c in_buffer, \c out_buffer, or \c in_out_buffers.
 * @tparam Arg Arguments given during kernel creation that specify the kernel
 * signature.
 * @tparam Args Arguments given during kernel creation that specify the kernel
 * signature.
 * @tparam CallArg First argument type used to call the kernel
 * @tparam CallArgs Other argument types used to call the kernel.
 * @param new_event The cl::Event generated involving the arguments.
 * @param m Arguments to the kernel that may be \c matrix_cls or not.
 * Non-matrices are ignored.
 * @param args Arguments to the kernel that may be matrices or not.
 * Non-matrices are ignored.
 */
template <typename T, require_same_t<T, cl::Event>* = nullptr>
inline void assign_events(const T&) {}

template <typename Arg, typename... Args, typename CallArg,
          typename... CallArgs>
inline void assign_events(const cl::Event& new_event, CallArg& m,
                          CallArgs&... args) {
  assign_event<Arg>(new_event, m);
  assign_events<Args...>(new_event, args...);
}

/** \ingroup kernel_executor_opencl
 * Select events from kernel arguments. Does nothing for non \c matrix_cl types.
 * @tparam T The argument type for a non \c matrix_cl, else the in/out/in_out
 * buffer types.
 * @tparam K The type of the \c matrix_cl
 * @param m If an \c matrix_cl, gets the event vector, else this argument does
 * nothing.
 * @return A vector of OpenCL events.
 */
template <typename T, require_not_matrix_cl_t<T>* = nullptr>
inline std::vector<cl::Event> select_events(const T& m) {
  return {};
}
template <typename T, typename K, require_matrix_cl_t<K>* = nullptr,
          require_same_t<T, in_buffer>* = nullptr>
inline const std::vector<cl::Event>& select_events(const K& m) {
  return m.write_events();
}
template <typename T, typename K, require_matrix_cl_t<K>* = nullptr,
          require_any_same_t<T, out_buffer, in_out_buffer>* = nullptr>
inline std::vector<cl::Event> select_events(K& m) {
  static_assert(!std::is_const<K>::value, "Can not write to const matrix_cl!");
  return m.read_write_events();
}

}  // namespace internal

/** \ingroup kernel_executor_opencl
 * Compile an OpenCL kernel.
 *
 * @param name The name for the kernel
 * @param sources A std::vector of strings containing the code for the kernel.
 * @param options The values of macros to be passed at compile time.
 */
inline auto compile_kernel(const char* name,
                           const std::vector<std::string>& sources,
                           const std::map<std::string, int>& options) {
  auto base_opts = opencl_context.base_opts();
  for (auto& it : options) {
    if (base_opts[it.first] > it.second) {
      base_opts[it.first] = it.second;
    }
  }
  std::string kernel_opts = "";
  for (auto&& comp_opts : base_opts) {
    kernel_opts += std::string(" -D") + comp_opts.first + "="
                   + std::to_string(comp_opts.second);
  }
  cl::Program program(opencl_context.context(), sources);
  try {
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

/** \ingroup kernel_executor_opencl
 * Creates functor for kernels
 *
 * @tparam Args Parameter pack of all kernel argument types.
 */
template <typename... Args>
struct kernel_cl {
 private:
  const char* name_;
  std::vector<std::string> sources_;
  std::map<std::string, int> opts_;
  mutable cl::Kernel kernel_;

 public:
  /** \ingroup kernel_executor_opencl
   * Creates functor for kernels that only need access to defining
   *  the global work size.
   * @param name The name for the kernel
   * @param sources A std::vector of strings containing the code for the kernel.
   * @param options The values of macros to be passed at compile time.
   */
  kernel_cl(const char* name, std::vector<std::string> sources,
            std::map<std::string, int> options = {})
      : name_(name), sources_(std::move(sources)), opts_(std::move(options)) {}

  /** \ingroup kernel_executor_opencl
   * Executes a kernel
   * @tparam CallArgs The types of the callee arguments.
   * @tparam Args Parameter pack of all kernel argument types.
   * @param global_thread_size The global work size.
   * @param args The arguments to pass to the kernel.
   * @return An Opencl event.
   */
  template <typename... CallArgs>
  auto operator()(cl::NDRange global_thread_size, CallArgs&&... args) const {
    if (kernel_() == NULL) {
      kernel_ = compile_kernel(name_, sources_, opts_);
      opencl_context.register_kernel_cache(&kernel_);
    }
    cl::EnqueueArgs eargs(opencl_context.queue(),
                          vec_concat(internal::select_events<Args>(args)...),
                          global_thread_size);
    cl::KernelFunctor<internal::to_const_buffer_t<Args>&...> kernel_functor(
        kernel_);
    cl::Event kern_event
        = kernel_functor(eargs, internal::get_kernel_args(args)...);
    internal::assign_events<Args...>(kern_event, args...);
    return kern_event;
  }

  /** \ingroup kernel_executor_opencl
   * Executes a kernel
   * @tparam CallArgs The types of the callee arguments.
   * @tparam Args Parameter pack of all kernel argument types.
   * @param global_thread_size The global work size.
   * @param thread_block_size The thread block size.
   * @param args The arguments to pass to the kernel.
   * @return An Opencl event.
   */
  template <typename... CallArgs>
  auto operator()(cl::NDRange global_thread_size, cl::NDRange thread_block_size,
                  CallArgs&&... args) const {
    if (kernel_() == NULL) {
      kernel_ = compile_kernel(name_, sources_, opts_);
      opencl_context.register_kernel_cache(&kernel_);
    }
    cl::EnqueueArgs eargs(opencl_context.queue(),
                          vec_concat(internal::select_events<Args>(args)...),
                          global_thread_size, thread_block_size);
    cl::KernelFunctor<internal::to_const_buffer_t<Args>&...> kernel_functor(
        kernel_);
    cl::Event kern_event
        = kernel_functor(eargs, internal::get_kernel_args(args)...);
    internal::assign_events<Args...>(kern_event, args...);
    return kern_event;
  }

  /** \ingroup kernel_executor_opencl
   * Retrieves an option used for compiling the kernel.
   * @param option_name which option to retrieve
   * @return option value
   */
  int get_option(const std::string option_name) const {
    return std::min(opts_.at(option_name),
                    opencl_context.base_opts().at(option_name));
  }
};

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
