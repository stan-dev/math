#ifndef STAN_MATH_REV_CORE_PRECOMPUTED_GRADIENTS_HPP
#define STAN_MATH_REV_CORE_PRECOMPUTED_GRADIENTS_HPP

#include <stan/math/prim/err/check_consistent_sizes.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>
// #include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/var.hpp>
#include <algorithm>
#include <vector>
#include <tuple>

namespace stan {
namespace math {

/**
 * A variable implementation taking a sequence of operands and
 * partial derivatives with respect to the operands.
 *
 * Stan users should use function precomputed_gradients()
 * directly.
 *
 * @tparam ContainerOperands tuple of any container operands (var_value
 * containing Eigen types)
 * @tparam ContainerGradients tupleof any container gradients (Eigen types)
 */
template <typename ContainerOperands = std::tuple<>,
          typename ContainerGradients = std::tuple<>>
class precomputed_gradients_vari_template : public vari {
 protected:
  const size_t size_;
  vari** varis_;
  double* gradients_;
  static_assert(std::tuple_size<ContainerOperands>::value
                    == std::tuple_size<ContainerGradients>::value,
                "precomputed_gradients_vari: ContainerOperands and "
                "ContainerGradients should have same size!");
  static constexpr size_t N_containers
      = std::tuple_size<ContainerOperands>::value;
  ContainerOperands container_operands_;
  ContainerGradients container_gradients_;

  /**
   * Checks that sizes of containers in container operands and
   * container_gradients match.
   *
   * @throws std::invalid_argument sizes do not match
   */
  template <size_t... Is>
  void check_sizes(std::index_sequence<Is...>) {
    static_cast<void>(std::initializer_list<int>{
        (check_size_match("precomputed_gradients_vari", "rows of operands",
                          std::get<Is>(container_operands_).vi_->rows(),
                          "rows of gradients",
                          std::get<Is>(container_gradients_).rows()),
         check_size_match("precomputed_gradients_vari", "cols of operands",
                          std::get<Is>(container_operands_).vi_->cols(),
                          "cols of gradients",
                          std::get<Is>(container_gradients_).cols()),
         0)...});
  }

 public:
  /**
   * Construct a precomputed vari with the specified value,
   * operands, gradients and optionally container operands and containers of
   * gradients.
   *
   * @tparam COps tuple of any container operands (var_value
   * containing Eigen types)
   * @tparam CGrads tupleof any container gradients (Eigen types)
   * @param[in] val The value of the variable.
   * @param[in] size Size of operands and gradients
   * @param[in] varis Operand implementations.
   * @param[in] gradients Gradients with respect to operands.
   * @param container_operands any container operands
   * @param container_gradients any container gradients
   */
  template <typename COps = std::tuple<>, typename CGrads = std::tuple<>>
  precomputed_gradients_vari_template(double val, size_t size, vari** varis,
                                      double* gradients,
                                      COps&& container_operands
                                      = std::tuple<>(),
                                      CGrads&& container_gradients
                                      = std::tuple<>())
      : vari(val),
        size_(size),
        varis_(varis),
        gradients_(gradients),
        container_operands_(std::forward<COps>(container_operands)),
        container_gradients_(std::forward<CGrads>(container_gradients)) {
    check_sizes(std::make_index_sequence<N_containers>());
  }

  /**
   * Construct a precomputed vari with the specified value,
   * operands, gradients and optionally container operands and containers of
   * gradients.
   *
   * @tparam Arith An arithmetic type
   * @tparam VecVar A vector of vars
   * @tparam VecArith A vector of arithmetic types
   * @tparam COps tuple of any container operands (var_value
   * containing Eigen types)
   * @tparam CGrads tupleof any container gradients (Eigen types)
   * @param[in] val The value of the variable.
   * @param[in] vars Vector of operands.
   * @param[in] gradients Vector of partial derivatives of value
   * with respect to operands.
   * @param container_operands any container operands
   * @param container_gradients any container gradients
   * @throws std::invalid_argument if the sizes of the vectors
   * don't match.
   */
  template <typename Arith, typename VecVar, typename VecArith,
            typename COps = std::tuple<>, typename CGrads = std::tuple<>,
            require_all_vector_t<VecVar, VecArith>* = nullptr>
  precomputed_gradients_vari_template(Arith val, const VecVar& vars,
                                      const VecArith& gradients,
                                      COps&& container_operands
                                      = std::tuple<>(),
                                      CGrads&& container_gradients
                                      = std::tuple<>())
      : vari(val),
        size_(vars.size()),
        varis_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            vars.size())),
        gradients_(ChainableStack::instance_->memalloc_.alloc_array<double>(
            vars.size())),
        container_operands_(std::forward<COps>(container_operands)),
        container_gradients_(std::forward<CGrads>(container_gradients)) {
    check_consistent_sizes("precomputed_gradients_vari", "vars", vars,
                           "gradients", gradients);
    check_sizes(std::make_index_sequence<N_containers>());
    for (size_t i = 0; i < vars.size(); ++i) {
      varis_[i] = vars[i].vi_;
    }
    std::copy(gradients.begin(), gradients.end(), gradients_);
  }

  /**
   * Implements the chain rule for this variable, using the
   * prestored operands and gradient.
   */
  void chain() {
    for (size_t i = 0; i < size_; ++i) {
      varis_[i]->adj_ += adj_ * gradients_[i];
    }
    index_apply<N_containers>([this](auto... Is) {
      static_cast<void>(std::initializer_list<int>{
          (std::get<Is>(this->container_operands_).adj()
           += this->adj_ * std::get<Is>(this->container_gradients_),
           0)...});
    });
  }
};

using precomputed_gradients_vari
    = precomputed_gradients_vari_template<std::tuple<>, std::tuple<>>;

/**
 * This function returns a var for an expression that has the
 * specified value, vector of operands, and vector of partial
 * derivatives of value with respect to the operands.
 *
 * @tparam Arith An arithmetic type
 * @tparam VecVar A vector of vars
 * @tparam VecArith A vector of arithmetic types
 * @tparam ContainerOperands tuple of any container operands (var_value
 * containing Eigen types)
 * @tparam ContainerGradients tupleof any container gradients (Eigen types)
 * @param[in] value The value of the resulting dependent variable.
 * @param[in] operands operands.
 * @param[in] gradients vector of partial derivatives of result with
 * respect to operands.
 * @param container_operands any container operands
 * @param container_gradients any container gradients
 * @return An autodiff variable that uses the precomputed
 * gradients provided.
 */
template <typename Arith, typename VecVar, typename VecArith,
          typename ContainerOperands = std::tuple<>,
          typename ContainerGradinets = std::tuple<>>
inline var precomputed_gradients(Arith value, const VecVar& operands,
                                 const VecArith& gradients,
                                 ContainerOperands&& container_operands
                                 = std::tuple<>(),
                                 ContainerGradinets&& container_gradients
                                 = std::tuple<>()) {
  return {
      new precomputed_gradients_vari_template<std::decay_t<ContainerOperands>,
                                              std::decay_t<ContainerGradinets>>(
          value, operands, gradients,
          std::forward<ContainerOperands>(container_operands),
          std::forward<ContainerGradinets>(container_gradients))};
}

}  // namespace math
}  // namespace stan
#endif
