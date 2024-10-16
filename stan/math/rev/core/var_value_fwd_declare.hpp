#ifndef STAN_MATH_REV_CORE_VAR_VALUE_FWD_DECLARE_HPP
#define STAN_MATH_REV_CORE_VAR_VALUE_FWD_DECLARE_HPP

#include <type_traits>

namespace stan {
namespace internal {
template <typename T, typename = void, typename = void>
struct arena_type_impl {
  static_assert(1, "This should never be instantiated");
  using type = void;
};
}

/**
 * Determines a type that can be used in place of `T` that does any dynamic
 * allocations on the AD stack. This way resulting types are trivially
 * destructible and can be used in vari classes.
 */
template <typename T>
using arena_t = typename internal::arena_type_impl<std::decay_t<T>>::type;

namespace math {



// forward declaration of var
// forward declare
template <typename Vari>
static void grad(Vari* vi);

template <typename T, typename = void>
class vari_value;
using vari = vari_value<double, void>;

template <typename T, typename = void>
class var_value;
using var = var_value<double, void>;

class chainable_alloc;

template <typename T>
struct arena_allocator;

/**
 * Equivalent to `Eigen::Matrix`, except that the data is stored on AD stack.
 * That makes these objects trivially destructible and usable in `vari`s.
 *
 * @tparam MatrixType Eigen matrix type this works as (`MatrixXd`, `VectorXd`,
 * ...)
 */
template <typename MatrixType, typename = void>
class arena_matrix;


/**
 * The type of a matrix holding <code>var</code>
 * values.
 */
using matrix_v = Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * The type of a (column) vector holding <code>var</code>
 * values.
 */
using vector_v = Eigen::Matrix<var, Eigen::Dynamic, 1>;

/**
 * The type of a row vector holding <code>var</code>
 * values.
 */
using row_vector_v = Eigen::Matrix<var, 1, Eigen::Dynamic>;

/**
 * The type of a matrix holding <code>vari*</code>
 * values.
 */
using matrix_vi = Eigen::Matrix<vari*, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * The type of a (column) vector holding <code>vari*</code>
 * values.
 */
using vector_vi = Eigen::Matrix<vari*, Eigen::Dynamic, 1>;

/**
 * The type of a row vector holding <code>vari*</code>
 * values.
 */
using row_vector_vi = Eigen::Matrix<vari*, 1, Eigen::Dynamic>;

}  // namespace math
}  // namespace stan
#endif
