#ifndef STAN_MATH_OPENCL_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_TRIANGULAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <Eigen/Core>
#include <type_traits>

namespace stan {
namespace math {
enum class matrix_cl_view { Diagonal = 0, Lower = 1, Upper = 2, Entire = 3 };

/**
 * Performs bitwise @c or to deduce return type from adding two @c matrix_cls
 * together
 * @param left_view first view
 * @param right_view second view
 * @return combined view
 */
inline const matrix_cl_view operator+(const matrix_cl_view& left_view,
                                      const matrix_cl_view& right_view) {
  typedef typename std::underlying_type<matrix_cl_view>::type underlying;
  return static_cast<matrix_cl_view>(static_cast<underlying>(left_view)
                                     | static_cast<underlying>(right_view));
}

/**
 * Performs bitwise @c and to deduce return type from adding two @c matrix_cls
 * together.
 * @param left_view first view
 * @param right_view second view
 * @return common nonzero part
 */
inline const matrix_cl_view operator*(const matrix_cl_view& left_view,
                                      const matrix_cl_view& right_view) {
  typedef typename std::underlying_type<matrix_cl_view>::type underlying;
  return static_cast<matrix_cl_view>(static_cast<underlying>(left_view)
                                     & static_cast<underlying>(right_view));
}

/**
 * Check whether a view contains certain nonzero part
 * @param left_view view to check
 * @param right_view part to check for (usually `Lower` or `Upper`)
 * @return true, if `a` has part `b` nonzero
 */
inline bool is_not_diagonal(const matrix_cl_view& left_view,
                            const matrix_cl_view& right_view) {
  return static_cast<bool>(left_view * right_view);
}

/**
 * Transposes a triangular view - swaps lower and upper parts.
 * @param view_type view to transpose
 * @return transposition of input
 */
inline const matrix_cl_view transpose(const matrix_cl_view& view_type) {
  if (view_type == matrix_cl_view::Lower) {
    return matrix_cl_view::Upper;
  }
  if (view_type == matrix_cl_view::Upper) {
    return matrix_cl_view::Lower;
  }
  return view_type;
}

/**
 * Inverts a triangular view. Parts that are zero in the input become nonzero in
 * output and vice versa.
 * @param view_type view to invert
 * @return inverted view
 */
inline const matrix_cl_view invert(const matrix_cl_view& view_type) {
  typedef typename std::underlying_type<matrix_cl_view>::type underlying;
  return static_cast<matrix_cl_view>(
      static_cast<underlying>(matrix_cl_view::Entire)
      & ~static_cast<underlying>(view_type));
}

/**
 * Creates a triangular view from `Eigen::UpLoType`. `Eigen::Lower`,
 * `Eigen::StrictlyLower` and `Eigen::UnitLower` become
 * `PartialViewCL::Lower`. Similar for `Upper`. Any other view becomes
 * `PartialViewCL::Entire`.
 * @param eigen_type `UpLoType` to create a view from
 * @return triangular view
 */
inline matrix_cl_view from_eigen_triangular_type(Eigen::UpLoType eigen_type) {
  if (eigen_type & Eigen::Lower) {
    return matrix_cl_view::Lower;
  }
  if (eigen_type & Eigen::Upper) {
    return matrix_cl_view::Upper;
  }
  return matrix_cl_view::Entire;
}

enum class TriangularMapCL { UpperToLower = 0, LowerToUpper = 1 };

// \cond
static const char* view_kernel_helpers = STRINGIFY(
    // \endcond
    /**
     * Combines two triangular views. Result is nonzero, where any of the inputs
     * is nonzero.
     * @param a first view
     * @param b second view
     * @return combined view
     */
    int combine(int a, int b) { return a | b; }

    /**
     * Determines common nonzero part of the inputs. Result is nonzero, where
     * both inputs are nonzero.
     * @param a first view
     * @param b second view
     * @return common nonzero part
     */
    int commonNonzeroPart(int a, int b) { return a & b; }

    /**
     * Check whether a view contains certain nonzero part
     * @param a view to check
     * @param b part to check for (usually `Lower` or `Upper`)
     * @return true, if `a` has part `b` nonzero
     */
    bool containsNonzeroPart(int a, int b) { return commonNonzeroPart(a, b); }
    // \cond
);
// \endcond

}  // namespace math
}  // namespace stan
#endif
#endif
