#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_holder_cl__HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_holder_cl__HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/wrapper.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <memory>
#include <string>
#include <type_traits>
#include <set>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents a no-op in kernel generator expressions. This object also owns
 * pointers to dynamically allocated objects used in its argument expression.
 * When this object is destructed, those objects are deleted.
 * @tparam Derived derived type
 * @tparam T type of the argument
 */
template <typename T, typename... Ptrs>
class holder_cl_
    : public operation_cl<holder_cl_<T, Ptrs...>,
                          typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<holder_cl_<T, Ptrs...>, Scalar, T>;
  using base::var_name;

 protected:
  std::tuple<std::unique_ptr<Ptrs>...> m_unique_ptrs;

 public:
  /**
   * Constructor
   * @param a argument expression
   * @param ptrs pointers this object takes ownership of
   */
  holder_cl_(T&& a, Ptrs*... ptrs)
      : base(std::forward<T>(a)),
        m_unique_ptrs(std::unique_ptr<Ptrs>(ptrs)...) {}

};

/**
 * Constructs a no-op operation that also holds pointer to some other
 * expressions, allocated on heap. When the object is destructed those
 * expressions will be deleted.
 * @tparam T type of argument expression
 * @tparam Ptrs types of pointers
 * @param a argument expression
 * @param ptrs pointers to objects the constructed object will own.
 * @return holder_cl_ operation
 */
template <typename T, typename... Ptrs>
auto holder_cl(T&& a, Ptrs*... ptrs) {
  return holder_cl_<as_operation_cl_t<T>, Ptrs...>(
      as_operation_cl(std::forward<T>(a)), ptrs...);
}

namespace internal {
template <typename T>
auto holder_cl_handle_element(const T& a, const T*& res) {
  res = &a;
  return std::make_tuple();
}

template <typename T>
auto holder_cl_handle_element(std::remove_reference_t<T>&& a, const T*& res) {
  res = new T(std::move(a));
  return std::make_tuple(res);
}

template <typename T, std::size_t... Is, typename... Args>
auto make_holder_cl_impl2(T&& expr, std::index_sequence<Is...>,
                          const std::tuple<Args*...>& ptrs) {
  return holder_cl(std::forward<T>(expr), std::get<Is>(ptrs)...);
}

template <typename T, std::size_t... Is, typename... Args>
auto make_holder_cl_impl(const T& func, std::index_sequence<Is...>,
                         Args&&... args) {
  std::tuple<const std::remove_reference_t<Args>*...> res;
  auto ptrs = std::tuple_cat(
      holder_cl_handle_element(std::forward<Args>(args), std::get<Is>(res))...);
  return make_holder_cl_impl2(
      func(*std::get<Is>(res)...),
      std::make_index_sequence<std::tuple_size<decltype(ptrs)>::value>(), ptrs);
}

}  // namespace internal

/**
 * Constructs an expression from given arguments using given functor.
 * This is similar to calling the functor with given arguments. Except that any
 * rvalue argument will be moved to heap first. The arguments moved to heap are
 * deleted once the expression is destructed.
 * @tparam T type of functor
 * @tparam Args types of arguments
 * @param func the functor
 * @param args arguments for the functor
 */
template <typename T, typename... Args>
auto make_holder_cl(const T& func, Args&&... args) {
  return internal::make_holder_cl_impl(
      func, std::make_index_sequence<sizeof...(Args)>(),
      std::forward<Args>(args)...);
}

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
