#ifndef STAN_MATH_PRIM_META_APPEND_RETURN_TYPE_HPP
#define STAN_MATH_PRIM_META_APPEND_RETURN_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <vector>

namespace stan {
namespace math {
/** \ingroup type_trait
 * This template metaprogram is used to compute the return type for
 * append_array
 *
 * This base implementation assumes the template parameters are not
 * vector types and calculates their return type using
 * <code>return_type</code>.
 *
 * @tparam T1 First type to be promoted
 * @tparam T2 Second type to be promoted
 */
template <typename T1, typename T2>
struct append_return_type {
  using type = return_type_t<T1, T2>;
};

/** \ingroup type_trait
 * This template metaprogram is used to compute the return type for
 * append_array
 *
 * If both types are ints, the type member is an int
 *
 * @tparam T1 First type to be promoted
 * @tparam T2 Second type to be promoted
 */
template <>
struct append_return_type<int, int> {
  using type = int;
};

/** \ingroup type_trait
 * This template metaprogram is used to compute the return type for
 * append_array
 *
 * If both types are Eigen::Matrices with the same Row/Column specification,
 * then the type member is another Eigen::Matrix with the same Row/Column
 * specification and a scalar type promoted using <code>return_type</code>.
 * Part of return type promotion logic for append_array
 *
 * @tparam T1 Scalar type of first matrix argument
 * @tparam T2 Scalar type of first matrix argument
 * @tparam R Eigen RowsAtCompileTime of both matrices
 * @tparam C Eigen ColsAtCompileTime of both matrices
 */
template <typename T1, typename T2, int R, int C>
struct append_return_type<Eigen::Matrix<T1, R, C>, Eigen::Matrix<T2, R, C> > {
  using type = typename Eigen::Matrix<return_type_t<T1, T2>, R, C>;
};

/** \ingroup type_trait
 * This template metaprogram is used to compute the return type for
 * append_array
 *
 * If the types of both template arguments are std::vectors, the type member
 * is recursively computed as the append_return_type of the scalar types
 * associated with those std::vectors.
 *
 * @tparam T1 Element type of first std::vector
 * @tparam T2 Element type of second std::vector
 */
template <typename T1, typename T2>
struct append_return_type<std::vector<T1>, std::vector<T2> > {
  using type = typename std::vector<typename append_return_type<T1, T2>::type>;
};
}  // namespace math
}  // namespace stan
#endif
