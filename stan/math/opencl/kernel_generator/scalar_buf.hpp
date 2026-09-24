#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_SCALAR_BUF_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_SCALAR_BUF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err/invalid_argument.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */
/**
 * Represents a device scalar (`opencl::ScalarCl<double>`) in kernel generator
 * expressions. Every thread reads the first element of the 1x1 buffer, so the
 * scalar broadcasts against expressions of any size.
 * @tparam T the 1x1 \c matrix_cl backing the device scalar, or a reference to
 * it
 */
template <typename T>
class scalar_buf_ : public operation_cl<scalar_buf_<T>, double> {
 protected:
  T a_;

  /**
   * Key identifying the backing buffer in the maps of already generated
   * operations. A `load_` of the same `matrix_cl` is keyed by the matrix's
   * address and generates different kernel arguments, so a different address
   * inside the same matrix is used here.
   * @return key for this buffer
   */
  inline const void* key() const noexcept { return &a_.read_events(); }

 public:
  using Scalar = double;
  using base = operation_cl<scalar_buf_<T>, Scalar>;
  using base::var_name_;
  static_assert(
      std::is_same<std::decay_t<T>, matrix_cl<double>>::value,
      "scalar_buf_: argument must be the matrix_cl<double> of a ScalarCl!");

  /**
   * Constructor
   * @param a 1x1 \c matrix_cl backing a device scalar
   */
  explicit scalar_buf_(T&& a) : a_(std::forward<T>(a)) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline scalar_buf_<T&> deep_copy() & { return scalar_buf_<T&>(a_); }
  inline scalar_buf_<const T&> deep_copy() const& {
    return scalar_buf_<const T&>(a_);
  }
  inline scalar_buf_<T> deep_copy() && {
    return scalar_buf_<T>(std::forward<T>(a_));
  }

  /**
   * Generates kernel code for this expression.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param name_gen name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether caller already handled matrix view
   * @return part of kernel with code for this expression
   */
  inline kernel_parts get_kernel_parts(
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name, bool view_handled) const {
    kernel_parts res{};
    const void* k = key();
    if (generated.count(k) == 0) {
      const char* arg_var_name;
      this->var_name_ = name_gen.generate();
      if (generated_all.count(k) == 0) {
        generated_all[k] = this->var_name_.c_str();
        arg_var_name = this->var_name_.c_str();
        res.args
            = "__global " + type_str<Scalar>() + "* " + var_name_ + "_global, ";
      } else {
        arg_var_name = generated_all[k];
      }
      generated[k] = this->var_name_.c_str();
      res.body = type_str<Scalar>() + " " + var_name_ + " = " + arg_var_name
                 + "_global[0];\n";
    } else {
      this->var_name_ = generated[k];
    }
    return res;
  }

  /**
   * Sets kernel arguments for this expression.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      cl::Kernel& kernel, int& arg_num) const {
    const void* k = key();
    if (generated_all.count(k) == 0) {
      generated_all[k] = "";
      kernel.setArg(arg_num++, a_.buffer());
    }
  }

  /**
   * Adds read event to the buffer used in this expression.
   * @param e the event to add
   */
  inline void add_read_event(cl::Event& e) const { a_.add_read_event(e); }

  /**
   * Adds all write events on the buffer used by this expression to a list.
   * @param[out] events List of all events.
   */
  inline void get_write_events(std::vector<cl::Event>& events) const {
    events.insert(events.end(), a_.write_events().begin(),
                  a_.write_events().end());
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression. A device scalar broadcasts, so this is dynamic.
   * @return number of rows
   */
  inline int rows() const { return base::dynamic; }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression. A device scalar broadcasts, so this is dynamic.
   * @return number of columns
   */
  inline int cols() const { return base::dynamic; }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
  }

  /**
   * Collects data that is needed beside types to uniquely identify a kernel
   * generator expression.
   * @param[out] uids ids of unique matrix accesses
   * @param[in,out] id_map map from memory addresses to unique ids
   * @param[in,out] next_id next unique id to use
   */
  inline void get_unique_matrix_accesses(
      std::vector<int>& uids, std::unordered_map<const void*, int>& id_map,
      int& next_id) const {
    const void* k = key();
    if (id_map.count(k) == 0) {
      id_map[k] = next_id;
      uids.push_back(next_id);
      next_id++;
    } else {
      uids.push_back(id_map[k]);
    }
  }
};

/**
 * Represents an expression whose result is a single element, used to evaluate
 * an expression into a device scalar. Scalar-only expressions have dynamic size
 * and are evaluated with a single thread. Expressions with a size other than
 * 1x1 are rejected.
 * @tparam T type of the argument
 */
template <typename T>
class scalar_result_
    : public operation_cl<scalar_result_<T>,
                          typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<scalar_result_<T>, Scalar, T>;
  using base::var_name_;

  /**
   * Constructor
   * @param a expression
   * @throw std::invalid_argument if the expression has a size other than 1x1
   */
  explicit scalar_result_(T&& a) : base(std::forward<T>(a)) {
    const auto& arg = this->template get_arg<0>();
    if (arg.rows() != base::dynamic && arg.rows() != 1) {
      invalid_argument("ScalarCl", "Rows of the assigned expression",
                       arg.rows(), "are ", ", but must be 1");
    }
    if (arg.cols() != base::dynamic && arg.cols() != 1) {
      invalid_argument("ScalarCl", "Columns of the assigned expression",
                       arg.cols(), "are ", ", but must be 1");
    }
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return scalar_result_<std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy)};
  }

  /**
   * Number of rows of the result.
   * @return 1
   */
  inline int rows() const { return 1; }

  /**
   * Number of columns of the result.
   * @return 1
   */
  inline int cols() const { return 1; }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
  }
};
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
