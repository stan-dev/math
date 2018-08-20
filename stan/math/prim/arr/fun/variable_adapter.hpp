#ifndef STAN_MATH_PRIM_ARR_FUN_VARIABLE_ADAPTER_HPP
#define STAN_MATH_PRIM_ARR_FUN_VARIABLE_ADAPTER_HPP

#include <stan/math/prim/scal/err/check_less.hpp>
#include <vector>
#include <cstddef>
#include <iostream>

namespace stan {
namespace math {

template <typename... Targs>
class variable_adapter {
 public:
  std::tuple<Targs...> args_;

  size_t size_;

 protected:
  template <typename T, typename... Pargs>
  size_t count_memory(size_t count, const std::vector<T>& x,
                      const Pargs&... args) {
    return count_memory(count + x.size(), args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const std::vector<int>& x,
                      const Pargs&... args) {
    return count_memory(count, args...);
  }

  template <typename T, typename... Pargs>
  size_t count_memory(size_t count, const T& x, const Pargs&... args) {
    return count_memory(count + 1, args...);
  }

  template <typename... Pargs>
  size_t count_memory(size_t count, const int& x, const Pargs&... args) {
    return count_memory(count, args...);
  }

  size_t count_memory(size_t count) { return count; }

  template <typename T>
  T& get(int i, std::vector<T>& arg) {
    return arg[i];
  }

  template <typename T>
  T& get(int i, T& arg) {
    return arg;
  }

  template <typename T, typename... Pargs>
  T& get(int i, std::vector<T>& arg, Pargs&... args) {
    if (i < arg.size())
      return arg[i];
    else
      return get(i - arg.size(), args...);
  }

  template <typename... Pargs>
  auto& get(int i, std::vector<int>& arg, Pargs&... args) {
    return get(i, args...);
  }

  template <typename T, typename... Pargs>
  T& get(int i, T& arg, Pargs&... args) {
    if (i == 0)
      return arg;
    else
      return get(i - 1, args...);
  }

  template <typename... Pargs>
  auto& get(int i, int& arg, Pargs&... args) {
    return get(i, args...);
  }

 public:
  variable_adapter(const Targs&... args)
      : args_(std::make_tuple(args...)), size_(count_memory(0, args...)) {}

  int size() const { return size_; }

  auto& operator()(int i) {
    check_less("variable_adapter::operator()", "i", i, size_);

    return std::apply(
        [ i, this ](auto&... args) -> auto& { return this->get(i, args...); },
        args_);
  }
};

}  // namespace math
}  // namespace stan
#endif
