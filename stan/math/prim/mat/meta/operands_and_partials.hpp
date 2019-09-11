#ifndef STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP
#define STAN_MATH_PRIM_MAT_META_OPERANDS_AND_PARTIALS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/broadcast_array.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/mat/meta/value_type.hpp>
#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {

/* This class will be used for both multivariate (nested container)
   operands_and_partials edges as well as for the univariate case.
 */
template <typename ViewElt, typename Op>
class ops_partials_edge<ViewElt, Op, std::enable_if_t<is_eigen<Op>::value && std::is_arithmetic<std::decay_t<value_type_t<Op>>>::value>> {
 public:
  using partials_t = empty_broadcast_array<ViewElt, Op>;
  partials_t partials_;
  empty_broadcast_array<partials_t, Op> partials_vec_;
  ops_partials_edge() {}
  explicit ops_partials_edge(const Op& /* ops */) {}

 private:
  template <typename, typename, typename, typename, typename, typename>
  friend class stan::math::operands_and_partials;

  void dump_partials(double* /* partials */) const {}  // reverse mode
  void dump_operands(void* /* operands */) const {}    // reverse mode
  double dx() const { return 0; }                      // used for fvars
  int size() const { return 0; }
};

template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<is_std_vector<Op>::value
                     && is_eigen<value_type_t<Op>>::value && std::is_arithmetic<std::decay_t<value_type_t<value_type_t<Op>>>>::value>> {
 public:
  using vector_scalar = value_type_t<Op>;
  using eigen_matrix = value_type_t<Op>;
  using partials_t = empty_broadcast_array<ViewElt, eigen_matrix>;
  empty_broadcast_array<partials_t, eigen_matrix> partials_vec_;
  ops_partials_edge() {}
  explicit ops_partials_edge(const Op& /* ops */) {}

 private:
  template <typename, typename, typename, typename, typename, typename,
            typename>
  friend class stan::math::operands_and_partials;

  void dump_partials(double* /* partials */) const {}  // reverse mode
  void dump_operands(void* /* operands */) const {}    // reverse mode
  double dx() const { return 0; }                      // used for fvars
  int size() const { return 0; }
};

template <typename ViewElt, typename Op>
class ops_partials_edge<
    ViewElt, Op,
    std::enable_if_t<is_std_vector<Op>::value
                     && is_std_vector<value_type_t<Op>>::value && std::is_arithmetic<std::decay_t<value_type_t<value_type_t<Op>>>>::value>> {
 public:
  using vector_scalar = value_type_t<Op>;
  using inner_vector_scalar = value_type_t<value_type_t<Op>>;
  using partials_t = empty_broadcast_array<ViewElt, Op>;
  partials_t partials_;
  empty_broadcast_array<partials_t, Op> partials_vec_;
  ops_partials_edge() {}
  explicit ops_partials_edge(const Op& /* ops */) {}

 private:
  template <typename, typename, typename, typename, typename, typename,
            typename>
  friend class stan::math::operands_and_partials;

  void dump_partials(double* /* partials */) const {}  // reverse mode
  void dump_operands(void* /* operands */) const {}    // reverse mode
  double dx() const { return 0; }                      // used for fvars
  int size() const { return 0; }
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
