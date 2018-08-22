#ifndef STAN_MATH_PRIM_ARR_FUN_VARIABLE_ADAPTER_HPP
#define STAN_MATH_PRIM_ARR_FUN_VARIABLE_ADAPTER_HPP

#include <stan/math/prim/scal/err/check_less.hpp>
#include <vector>
#include <cstddef>
#include <iostream>

namespace stan {
namespace math {

template <typename T, typename... Targs>
class variable_adapter {
 public:
  std::tuple<Targs...> args_;

  size_t size_;

 protected:
  template <typename... Pargs>
  size_t count_memory(size_t count, const std::vector<T>& x,
                      const Pargs&... args) {
    return count_memory(count + x.size(), args...);
  }

  template <typename R, typename... Pargs>
  size_t count_memory(size_t count, const std::vector<R>& x,
                      const Pargs&... args) {
    return count_memory(count, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const T& x, const Pargs&... args) {
    return count_memory(count + 1, args...);
  }

  template <typename R, typename... Pargs>
  size_t count_memory(size_t count, const R& x, const Pargs&... args) {
    return count_memory(count, args...);
  }

  size_t count_memory(size_t count) { return count; }

  T& get(size_t i, std::vector<T>& arg) { return arg[i]; }

  template <typename R>
  T& get(size_t i, std::vector<R>& arg) {
    throw std::runtime_error("This should not be reachable");
  }

  T& get(size_t i, T& arg) { return arg; }

  template <typename R>
  T& get(size_t i, R& arg) {
    throw std::runtime_error("This should not be reachable");
  }

  template <typename... Pargs>
  T& get(size_t i, std::vector<T>& arg, Pargs&... args) {
    if (i < arg.size())
      return arg[i];
    else
      return get(i - arg.size(), args...);
  }

  template <typename R, typename... Pargs>
  auto& get(size_t i, std::vector<R>& arg, Pargs&... args) {
    return get(i, args...);
  }

  template <typename... Pargs>
  T& get(size_t i, T& arg, Pargs&... args) {
    if (i == 0)
      return arg;
    else
      return get(i - 1, args...);
  }

  template <typename R, typename... Pargs>
  auto& get(size_t i, R& arg, Pargs&... args) {
    return get(i, args...);
  }

 public:
  variable_adapter(const Targs&... args)
      : args_(std::make_tuple(args...)), size_(count_memory(0, args...)) {}

  size_t size() const { return size_; }

  auto& operator()(size_t i) {
    check_less("variable_adapter::operator()", "i", i, size_);

    return std::apply(
        [ i, this ](auto&... args) -> auto& { return this->get(i, args...); },
        args_);
  }
};

template <typename T, typename... Targs>
auto variable_adapter_factory(const Targs&... args) {
  return variable_adapter<T, Targs...>(args...);
}

}  // namespace math
}  // namespace stan
#endif
