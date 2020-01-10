#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_OPERATION_LHS_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_OPERATION_LHS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <string>
#include <set>
#include <array>
#include <numeric>

namespace stan {
namespace math {

/**
 * Base for all kernel generator operations that can be used on left hand side
 * of an expression.
 * @tparam Derived derived type
 * @tparam Scalar scalar type of the result
 * @tparam Args types of arguments to this operation
 */
template <typename Derived, typename Scalar, typename... Args>
class operation_cl_lhs : public operation_cl<Derived, Scalar, Args...> {
 protected:
  using base = operation_cl<Derived, Scalar, Args...>;
  static constexpr int N = sizeof...(Args);
  using base::arguments_;
  using base::derived;

 public:
  using base::operation_cl;

  /**
   * generates kernel code for this expression if it appears on the left hand
   * side of an assigment.
   * @param[in,out] generated set of (pointer to) already generated operations
   * @param name_gen name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this expressions
   */
  inline kernel_parts get_kernel_parts_lhs(
      std::set<const operation_cl_base*>& generated, name_generator& name_gen,
      const std::string& i, const std::string& j) const {
    if (generated.count(this) == 0) {
      this->var_name = name_gen.generate();
    }
    std::string i_arg = i;
    std::string j_arg = j;
    this->derived().modify_argument_indices(i_arg, j_arg);
    std::array<kernel_parts, N> args_parts = index_apply<N>([&](auto... Is) {
      return std::array<kernel_parts, N>{
          std::get<Is>(this->arguments_)
              .get_kernel_parts_lhs(generated, name_gen, i_arg, j_arg)...};
    });
    kernel_parts res{};
    res.body = std::accumulate(
        args_parts.begin(), args_parts.end(), std::string(),
        [](const std::string& a, const kernel_parts& b) { return a + b.body; });
    if (generated.count(this) == 0) {
      generated.insert(this);
      res.args
          = std::accumulate(args_parts.begin(), args_parts.end(), std::string(),
                            [](const std::string& a, const kernel_parts& b) {
                              return a + b.args;
                            });
      kernel_parts my_part = index_apply<N>([&](auto... Is) {
        return this->derived().generate_lhs(
            i, j, std::get<Is>(this->arguments_).var_name...);
      });
      res.body += my_part.body;
      res.args += my_part.args;
    }
    return res;
  }

  /**
   * Sets view of the underlying matrix depending on which part is written.
   * @param bottom_diagonal Index of the top sub- or super- diagonal written
   * with nonzero elements.
   * @param top_diagonal Index of the top sub- or super- diagonal written with
   * nonzero elements.
   * @param bottom_zero_diagonal Index of the top sub- or super- diagonal
   * written with zeros if it ie more extreme than \c bottom_diagonal. Otherwise
   * it should be set to equal value as \c bottom_diagonal.
   * @param top_zero_diagonal Index of the top sub- or super- diagonal written
   * with zeros if it ie more extreme than \c top_diagonal. Otherwise it should
   * be set to equal value as \c top_diagonal.
   */
  inline void set_view(int bottom_diagonal, int top_diagonal,
                       int bottom_zero_diagonal, int top_zero_diagonal) const {
    index_apply<N>([&](auto... Is) {
      (void)std::initializer_list<int>{
          (std::get<Is>(this->arguments_)
               .set_view(bottom_diagonal, top_diagonal, bottom_zero_diagonal,
                         top_zero_diagonal),
           0)...};
    });
  }

  /**
   * Adds write event to any matrices used by nested expressions.
   * @param e the event to add
   */
  inline void add_write_event(cl::Event& e) const {
    index_apply<N>([&](auto... Is) {
      (void)std::initializer_list<int>{
          (std::get<Is>(this->arguments_).add_write_event(e), 0)...};
    });
  }
};

}  // namespace math
}  // namespace stan

#endif
#endif
