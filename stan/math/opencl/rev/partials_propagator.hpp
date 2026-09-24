#ifndef STAN_MATH_OPENCL_REV_PARTIALS_PROPAGATOR_HPP
#define STAN_MATH_OPENCL_REV_PARTIALS_PROPAGATOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/partials_propagator.hpp>
#include <stan/math/opencl/rev/arena_matrix_cl.hpp>
#include <stan/math/opencl/rev/scalar_cl.hpp>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Edge of an OpenCL partials propagator. Operands that are not autodiff have
 * no partials and receive nothing in the reverse pass.
 * @tparam Op type of the operand
 */
template <typename Op, typename Enable = void>
class ops_partials_edge_cl {
 public:
  template <typename T>
  explicit ops_partials_edge_cl(const T& /* op */) noexcept {}
  template <typename Adj>
  inline void propagate(const Adj& /* res_adj */) noexcept {}
};

/**
 * Edge for a CPU var operand. The partial is kept on the device. In the
 * reverse pass the product of the result adjoint and the partial is read back
 * and added to the CPU adjoint; this is the one transfer the CPU operand
 * requires.
 */
template <>
class ops_partials_edge_cl<var> {
 public:
  opencl::internal::scalar_cl_partial partials_;
  var operand_;
  explicit ops_partials_edge_cl(const var& op) : operand_(op) {}
  template <typename Adj>
  inline void propagate(const Adj& res_adj) {
    operand_.adj() += opencl::to_host(res_adj * partials_.value_);
  }
};

/**
 * Edge for a device var operand. The partial and the adjoint update stay on
 * the device.
 */
template <>
class ops_partials_edge_cl<opencl::ScalarCl<var>> {
 public:
  opencl::internal::scalar_cl_partial partials_;
  opencl::ScalarCl<var> operand_;
  explicit ops_partials_edge_cl(const opencl::ScalarCl<var>& op)
      : operand_(op) {}
  template <typename Adj>
  inline void propagate(const Adj& res_adj) {
    operand_.adj() += elt_multiply(as_operation_cl(res_adj),
                                   as_operation_cl(partials_.value_));
  }
};

/**
 * Edge for a var containing an OpenCL matrix. The partials are a matrix on the
 * device, and the adjoint update is a device expression broadcasting the
 * result adjoint.
 */
template <typename Op>
class ops_partials_edge_cl<var_value<Op>, require_kernel_expression_lhs_t<Op>> {
 public:
  arena_matrix_cl<value_type_t<Op>> partials_;
  var_value<Op> operand_;
  explicit ops_partials_edge_cl(const var_value<Op>& op)
      : partials_(constant(0.0, op.vi_->rows(), op.vi_->cols())),
        operand_(op) {}
  template <typename Adj>
  inline void propagate(const Adj& res_adj) {
    operand_.adj() += res_adj * partials_;
  }
};

/**
 * Partials propagator for OpenCL functions with at least one autodiff
 * operand. The result is a device var and the whole reverse pass runs on the
 * device, except for CPU var operands (see `ops_partials_edge_cl<var>`).
 */
template <typename ReturnType, typename... Ops>
class partials_propagator<ReturnType, require_rev_scalar_cl_t<ReturnType>,
                          Ops...> {
 public:
  std::tuple<ops_partials_edge_cl<plain_type_t<std::decay_t<Ops>>>...> edges_;

  template <typename... Types>
  explicit partials_propagator(Types&&... ops)
      : edges_(ops_partials_edge_cl<plain_type_t<std::decay_t<Ops>>>(
          std::forward<Types>(ops))...) {}

  /**
   * Builds the device var holding the value of the function. In the reverse
   * pass the adjoint of the result, times each operand's partials, is added
   * to the operand's adjoint.
   * @param value the value of the function
   * @return device var
   */
  inline opencl::ScalarCl<var> build(opencl::ScalarCl<double>&& value) {
    return opencl::make_callback_scalar_cl(
        std::move(value),
        [edges = std::move(edges_)](const auto& res_adj, const auto&) mutable {
          stan::math::for_each(
              [&res_adj](auto& edge) { edge.propagate(res_adj); }, edges);
        });
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
#endif
