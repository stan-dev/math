#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_MULTI_RESULT_KERNEL_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_MULTI_RESULT_KERNEL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/calc_if.hpp>
#include <stan/math/opencl/kernel_generator/check_cl.hpp>
#include <stan/math/opencl/kernel_generator/colwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/reduction_2d.hpp>
#include <stan/math/opencl/kernel_generator/load.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <algorithm>
#include <string>
#include <tuple>
#include <utility>
#include <map>
#include <vector>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */
namespace internal {

// Template parameter pack can only be at the end of the template list in
// structs. We need 2 packs for expressions and results, so we nest structs.
template <int N, typename... T_results>
struct multi_result_kernel_internal {
  template <typename... T_expressions>
  struct inner {
    static std::map<std::vector<int>, cl::Kernel> kernel_cache_;
    using next = typename multi_result_kernel_internal<
        N - 1, T_results...>::template inner<T_expressions...>;
    using T_current_result = std::remove_reference_t<
        std::tuple_element_t<N, std::tuple<T_results...>>>;
    using T_current_expression = std::remove_reference_t<
        std::tuple_element_t<N, std::tuple<T_expressions...>>>;
    /**
     * Generates list of all events kernel assigning expressions to results must
     * wait on. Also clears those events from matrices.
     * @param[out] events list of events
     * @param assignment_pairs pairs of result and expression
     */
    static void get_clear_events(
        std::vector<cl::Event>& events,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {
      next::get_clear_events(events, assignment_pairs);
      std::get<N>(assignment_pairs).second.get_write_events(events);
      std::get<N>(assignment_pairs).first.get_clear_read_write_events(events);
    }
    /**
     * Assigns the dimensions of expressions to matching results if possible.
     * Otherwise checks that dimensions match. Also checks that all expressions
     * require same number of threads.
     * @param n_rows number of threads in rows dimension of the first
     * expression
     * @param n_cols number of threads in rows dimension of the first
     * expression
     * @param assignment_pairs pairs of result and expression
     */
    static void check_assign_dimensions(
        int n_rows, int n_cols,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {
      next::check_assign_dimensions(n_rows, n_cols, assignment_pairs);
      const auto& expression = std::get<N>(assignment_pairs).second;
      const auto& result = std::get<N>(assignment_pairs).first;
      const char* function = "results.operator=";

      int expression_rows = expression.rows();
      int expression_cols = expression.cols();
      if (is_colwise_reduction<T_current_expression>::value
          && expression_cols == -1) {
        expression_cols = n_cols;
        expression_rows = expression.thread_rows();
        expression_rows = internal::colwise_reduction_wgs_rows(
            expression_rows < 0 ? n_rows : expression_rows, expression_cols);
      } else if (is_reduction_2d<T_current_expression>::value
                 && expression_cols == -1) {
        expression_rows = internal::colwise_reduction_wgs_rows(n_rows, n_cols);
        if (expression_rows == 0) {
          expression_cols = 0;
        } else {
          expression_cols = (n_cols + expression_rows - 1) / expression_rows;
        }
      } else {
        if (expression.thread_rows() != -1) {
          check_size_match(function, "Rows of ", "expression",
                           expression.thread_rows(), "rows of ",
                           "first expression", n_rows);
        } else {
          expression_rows = n_rows;
        }
        if (expression.thread_cols() != -1) {
          check_size_match(function, "Columns of ", "expression",
                           expression.thread_cols(), "columns of ",
                           "first expression", n_cols);
        } else {
          expression_cols = n_cols;
        }
      }
      result.check_assign_dimensions(expression_rows, expression_cols);
      int bottom_written = 1 - expression.rows();
      int top_written = expression.cols() - 1;
      std::pair<int, int> extreme_diagonals = expression.extreme_diagonals();
      result.set_view(std::max(extreme_diagonals.first, bottom_written),
                      std::min(extreme_diagonals.second, top_written),
                      bottom_written, top_written);
    }

    /**
     * Generates kernel source for assignment of expressions to results.
     * @param[in,out] generated map from (pointer to) already generated local
     * operations to variable names
     * @param[in,out] generated_all map from (pointer to) already generated all
     * operations to variable names
     * @param ng name generator
     * @param row_index_name variable name of the row index
     * @param col_index_name variable name of the column index
     * @param assignment_pairs pairs of result and expression
     * @return kernel parts for the kernel
     */
    static kernel_parts generate(
        std::map<const void*, const char*>& generated,
        std::map<const void*, const char*>& generated_all, name_generator& ng,
        const std::string& row_index_name, const std::string& col_index_name,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {
      kernel_parts parts
          = next::generate(generated, generated_all, ng, row_index_name,
                           col_index_name, assignment_pairs);

      kernel_parts parts0
          = std::get<N>(assignment_pairs)
                .second.get_whole_kernel_parts(
                    generated, generated_all, ng, row_index_name,
                    col_index_name, std::get<N>(assignment_pairs).first);
      parts += parts0;
      return parts;
    }

    /**
     * Sets kernel arguments.
     * @param[in,out] generated map from (pointer to) already generated local
     * operations to variable names
     * @param[in,out] generated_all map from (pointer to) already generated all
     * operations to variable names
     * @param kernel kernel to set arguments to
     * @param arg_num number of the next argument to set
     * @param assignment_pairs pairs of result and expression
     */
    static void set_args(
        std::map<const void*, const char*>& generated,
        std::map<const void*, const char*>& generated_all, cl::Kernel& kernel,
        int& arg_num,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {
      next::set_args(generated, generated_all, kernel, arg_num,
                     assignment_pairs);

      std::get<N>(assignment_pairs)
          .second.set_args(generated, generated_all, kernel, arg_num);
      std::get<N>(assignment_pairs)
          .first.set_args(generated, generated_all, kernel, arg_num);
    }

    /**
     * Adds event to matrices used in kernel.
     * @param e event to add
     * @param assignment_pairs pairs of result and expression
     */
    static void add_event(
        cl::Event e, const std::tuple<std::pair<T_results, T_expressions>...>&
                         assignment_pairs) {
      next::add_event(e, assignment_pairs);
      std::get<N>(assignment_pairs).second.add_read_event(e);
      std::get<N>(assignment_pairs).first.add_write_event(e);
    }

    /**
     * Collects data that is needed beside types to uniqly identify a kernel.
     * @param[out] uids ids of unique matrix accesses
     * @param[in,out] id_map map from memory addresses to unique ids
     * @param[in,out] next_id neqt unique id to use
     * @param assignment_pairs pairs of result and expression
     */
    static void get_unique_matrix_accesses(
        std::vector<int>& uids, std::map<const void*, int>& id_map,
        int& next_id,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {
      std::get<N>(assignment_pairs)
          .second.get_unique_matrix_accesses(uids, id_map, next_id);
      std::get<N>(assignment_pairs)
          .first.get_unique_matrix_accesses(uids, id_map, next_id);
      next::get_unique_matrix_accesses(uids, id_map, next_id, assignment_pairs);
    }
  };
};

// Specialization for N == -1 ends the recursion.
template <typename... T_results>
struct multi_result_kernel_internal<-1, T_results...> {
  template <typename... T_expressions>
  struct inner {
    static void get_clear_events(
        std::vector<cl::Event>& events,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {}

    static void check_assign_dimensions(
        int n_rows, int n_cols,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {}

    static kernel_parts generate(
        std::map<const void*, const char*>& generated,
        std::map<const void*, const char*>& generated_all, name_generator& ng,
        const std::string& row_index_name, const std::string& col_index_name,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {
      return {};
    }

    static void set_args(
        std::map<const void*, const char*>& generated,
        std::map<const void*, const char*>& generated_all, cl::Kernel& kernel,
        int& arg_num,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {}

    static void add_event(
        cl::Event e, const std::tuple<std::pair<T_results, T_expressions>...>&
                         assignment_pairs) {}

    static void get_unique_matrix_accesses(
        std::vector<int>& uids, std::map<const void*, int>& id_map,
        int& next_id,
        const std::tuple<std::pair<T_results, T_expressions>...>&
            assignment_pairs) {}
  };
};

template <int N, typename... T_results>
template <typename... T_expressions>
std::map<std::vector<int>, cl::Kernel> multi_result_kernel_internal<
    N, T_results...>::inner<T_expressions...>::kernel_cache_;

}  // namespace internal

/**
 * Represents multiple expressions that will be calculated in same kernel.
 * @tparam T_expressions types of expressions
 */
template <typename... T_expressions>
class expressions_cl {
 public:
  /**
   * Constructor.
   * @param expressions expressions that will be calculated in same kernel.
   */
  explicit expressions_cl(T_expressions&&... expressions)
      : expressions_(
          T_expressions(std::forward<T_expressions>(expressions))...) {}

 private:
  std::tuple<T_expressions...> expressions_;
  template <typename...>
  friend class results_cl;
  template <typename...>
  friend class adjoint_results_cl;
};

/**
 * Deduces types for constructing \c expressions_cl object.
 * @tparam T_expressions types of expressions
 * @param expressions expressions that will be used in same kernel.
 */
template <typename... T_expressions>
expressions_cl<T_expressions...> expressions(T_expressions&&... expressions) {
  return expressions_cl<T_expressions...>(
      std::forward<T_expressions>(expressions)...);
}

/**
 * Represents results that will be calculated in same kernel.
 * @tparam T_results types of results
 */
template <typename... T_results>
class results_cl {
 protected:
  std::tuple<T_results...> results_;

 public:
  /**
   * Constructor.
   * @param results results that will be calculated in same kernel
   */
  explicit results_cl(T_results&&... results)
      : results_(std::forward<T_results>(results)...) {}

  /**
   * Assigning \c expressions_cl object to \c results_ object generates and
   * executes the kernel that evaluates expressions and stores them into result
   * expressions this object contains.
   * @tparam T_expressions types of expressions
   * @param exprs expressions
   */
  template <typename... T_expressions,
            typename = std::enable_if_t<sizeof...(T_results)
                                        == sizeof...(T_expressions)>>
  void operator=(const expressions_cl<T_expressions...>& exprs) {
    index_apply<sizeof...(T_expressions)>([this, &exprs](auto... Is) {
      assignment_impl(std::tuple_cat(make_assignment_pair(
          std::get<Is>(results_), std::get<Is>(exprs.expressions_))...));
    });
  }

  /**
   * Incrementing \c results_ object by \c expressions_cl object
   * executes the kernel that evaluates expressions and increments results by
   * those expressions.
   * @tparam T_expressions types of expressions
   * @param exprs expressions
   */
  template <typename... T_expressions,
            typename = std::enable_if_t<sizeof...(T_results)
                                        == sizeof...(T_expressions)>>
  void operator+=(const expressions_cl<T_expressions...>& exprs) {
    index_apply<sizeof...(T_expressions)>([this, &exprs](auto... Is) {
      auto tmp = std::tuple_cat(make_assignment_pair(
          std::get<Is>(results_), std::get<Is>(exprs.expressions_))...);
      index_apply<std::tuple_size<decltype(tmp)>::value>(
          [this, &tmp](auto... Is2) {
            assignment_impl(std::make_tuple(std::make_pair(
                std::get<Is2>(tmp).first,
                std::get<Is2>(tmp).first + std::get<Is2>(tmp).second)...));
          });
    });
  }

  /**
   * Generates kernel source for evaluating given expressions into results held
   * by \c this.
   * @tparam T_expressions types of expressions
   * @param exprs expressions to generate kernel source for
   * @return kernel source
   */
  template <typename... T_expressions,
            typename = std::enable_if_t<sizeof...(T_results)
                                        == sizeof...(T_expressions)>>
  std::string get_kernel_source_for_evaluating(
      const expressions_cl<T_expressions...>& exprs) {
    return index_apply<sizeof...(T_expressions)>([this, &exprs](auto... Is) {
      return get_kernel_source_impl(std::tuple_cat(make_assignment_pair(
          std::get<Is>(results_), std::get<Is>(exprs.expressions_))...));
    });
  }

 protected:
  /**
   * Implementation of kernel source generation.
   * @tparam T_res types of results
   * @tparam T_expressions types of expressions
   * @param assignment_pairs pairs if result and expression
   */
  template <typename... T_res, typename... T_expressions>
  static std::string get_kernel_source_impl(
      const std::tuple<std::pair<T_res, T_expressions>...>& assignment_pairs) {
    using impl = typename internal::multi_result_kernel_internal<
        std::tuple_size<std::tuple<T_expressions...>>::value - 1,
        T_res...>::template inner<T_expressions...>;
    static const bool require_specific_local_size = std::max(
        {std::decay_t<T_expressions>::Deriv::require_specific_local_size...});

    name_generator ng;
    std::map<const void*, const char*> generated;
    std::map<const void*, const char*> generated_all;
    kernel_parts parts = impl::generate(generated, generated_all, ng, "i", "j",
                                        assignment_pairs);
    std::string src;
    if (require_specific_local_size) {
      src =
          parts.includes +
          "kernel void calculate(" + parts.args +
          "const int rows, const int cols){\n"
          "const int gid_i = get_global_id(0);\n"
          "const int gid_j = get_global_id(1);\n"
          "const int lid_i = get_local_id(0);\n"
          "const int lsize_i = get_local_size(0);\n"
          "const int gsize_i = get_global_size(0);\n"
          "const int gsize_j = get_global_size(1);\n"
          "const int wg_id_i = get_group_id(0);\n"
          "const int wg_id_j = get_group_id(1);\n"
          "const int n_groups_i = get_num_groups(0);\n"
          + parts.declarations +
          "for(int j = gid_j; j < cols; j+=gsize_j){\n"
          + parts.initialization +
          "for(int i = gid_i; i < rows; i+=gsize_i){\n"
          + parts.body
          + parts.body_suffix +
          "}\n"
          + parts.reduction_1d +
          "}\n"
          + parts.reduction_2d +
          "}\n";
    } else {
      src =
          parts.includes +
          "kernel void calculate(" +
          parts.args.substr(0, parts.args.size() - 2) +
          "){\n"
          "int i = get_global_id(0);\n"
          "int j = get_global_id(1);\n"
          + parts.declarations
          + parts.initialization
          + parts.body
          + parts.body_suffix
          + parts.reduction_1d
          + parts.reduction_2d +
          "}\n";
    }
    return src;
  }

  /**
   * Implementation of assignments of expressions to results
   * @tparam T_res types of results
   * @tparam T_expressions types of expressions
   * @param assignment_pairs pairs if result and expression
   */
  template <typename... T_res, typename... T_expressions>
  static void assignment_impl(
      const std::tuple<std::pair<T_res, T_expressions>...>& assignment_pairs) {
    using impl = typename internal::multi_result_kernel_internal<
        std::tuple_size<std::tuple<T_expressions...>>::value - 1,
        T_res...>::template inner<T_expressions...>;

    static const bool any_output = std::max(
        {false, !is_without_output<std::decay_t<T_expressions>>::value...});
    if (!any_output) {
      return;
    }

    static const bool require_specific_local_size = std::max(
        {std::decay_t<T_expressions>::Deriv::require_specific_local_size...});

    int n_rows = std::get<0>(assignment_pairs).second.thread_rows();
    int n_cols = std::get<0>(assignment_pairs).second.thread_cols();
    const char* function = "results_cl.assignment";
    impl::check_assign_dimensions(n_rows, n_cols, assignment_pairs);
    if (n_rows * n_cols == 0) {
      return;
    }
    if (n_rows < 0) {
      invalid_argument(function, "Number of rows of expression", n_rows,
                       " must be nonnegative, but is ",
                       " (broadcasted expressions can not be evaluated)");
    }
    if (n_cols < 0) {
      invalid_argument(function, "Number of columns of expression", n_cols,
                       " must be nonnegative, but is ",
                       " (broadcasted expressions can not be evaluated)");
    }

    std::vector<int> uids;
    std::map<const void*, int> id_map;
    int next_id = 0;
    impl::get_unique_matrix_accesses(uids, id_map, next_id, assignment_pairs);

    try {
      if (impl::kernel_cache_[uids]() == NULL) {
        std::string src = get_kernel_source_impl(assignment_pairs);
        auto opts = opencl_context.base_opts();
        opts["LOCAL_SIZE_"] = std::min(64, opts.at("LOCAL_SIZE_"));
        impl::kernel_cache_[uids] = opencl_kernels::compile_kernel(
            "calculate", {view_kernel_helpers, src}, opts);
        opencl_context.register_kernel_cache(&impl::kernel_cache_[uids]);
      }
      cl::Kernel& kernel = impl::kernel_cache_[uids];
      int arg_num = 0;

      std::map<const void*, const char*> generated;
      std::map<const void*, const char*> generated_all;
      impl::set_args(generated, generated_all, kernel, arg_num,
                     assignment_pairs);

      std::vector<cl::Event> events;
      impl::get_clear_events(events, assignment_pairs);
      cl::Event e;
      if (require_specific_local_size) {
        kernel.setArg(arg_num++, n_rows);
        kernel.setArg(arg_num++, n_cols);

        int local = std::min(64, opencl_context.base_opts().at("LOCAL_SIZE_"));

        int wgs_rows = internal::colwise_reduction_wgs_rows(n_rows, n_cols);
        int wgs_cols = (n_cols + wgs_rows - 1) / wgs_rows;

        opencl_context.queue().enqueueNDRangeKernel(
            kernel, cl::NullRange, cl::NDRange(local * wgs_rows, wgs_cols),
            cl::NDRange(local, 1), &events, &e);
      } else {
        opencl_context.queue().enqueueNDRangeKernel(kernel, cl::NullRange,
                                                    cl::NDRange(n_rows, n_cols),
                                                    cl::NullRange, &events, &e);
      }
      impl::add_event(e, assignment_pairs);
    } catch (const cl::Error& e) {
      check_opencl_error(function, e);
    }
  }

  /**
   * Implementation of assignments of no expressions to no results
   */
  static void assignment_impl(const std::tuple<>& /*assignment_pairs*/) {}

  /**
   * Makes a std::pair of one result and one expression and wraps it into a
   * tuple.
   * @param result result
   * @param expression expression
   * @return a tuple of pair of result and expression
   */
  template <typename T_result, typename T_expression,
            require_all_not_t<is_without_output<T_expression>,
                              conjunction<internal::is_scalar_check<T_result>,
                                          std::is_arithmetic<std::decay_t<
                                              T_expression>>>>* = nullptr>
  static auto make_assignment_pair(T_result&& result,
                                   T_expression&& expression) {
    return std::make_tuple(
        std::pair<as_operation_cl_t<T_result>, as_operation_cl_t<T_expression>>(
            as_operation_cl(std::forward<T_result>(result)),
            as_operation_cl(std::forward<T_expression>(expression))));
  }

  /**
   * If an expression does not need to be calculated this returns an empty tuple
   * @param result result
   * @param expression expression
   * @return a tuple of pair of result and expression
   */
  template <typename T_result, typename T_expression,
            require_t<is_without_output<T_expression>>* = nullptr>
  static auto make_assignment_pair(T_result&& result,
                                   T_expression&& expression) {
    return std::make_tuple();
  }

  /**
   * Checks on scalars are done separately in this overload instead of in
   * kernel.
   * @param result result - check
   * @param pass bool scalar
   * @return an empty tuple
   */
  template <typename T_check, typename T_pass,
            require_t<internal::is_scalar_check<T_check>>* = nullptr,
            require_integral_t<T_pass>* = nullptr>
  static std::tuple<> make_assignment_pair(T_check&& result, T_pass&& pass) {
    if (!pass) {
      std::stringstream s;
      s << result.function_ << ": " << result.err_variable_ << " = "
        << result.arg_.a_ << ", but it must be " << result.must_be_ << "!";
      throw std::domain_error(s.str());
    }
    return std::make_tuple();
  }
};

/**
 * Deduces types for constructing \c results_cl object.
 * @tparam T_results types of results
 * @param results results that will be calculated in same kernel.
 */
template <typename... T_results>
results_cl<T_results...> results(T_results&&... results) {
  return results_cl<T_results...>(std::forward<T_results>(results)...);
}

// an implementation that needs results and expressions defined (not just
// declared)
template <typename T>
void check_cl_<T>::operator=(bool condition) {
  results(*this) = expressions(condition);
}

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
