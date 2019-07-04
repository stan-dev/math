#ifndef STAN_MATH_OPENCL_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_TRIANGULAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <Eigen/Core>
#include <type_traits>

namespace stan {
namespace math {
enum class TriangularViewCL { Diagonal = 0, Lower = 1, Upper = 2, Entire = 3 };

/**
 * Combines two triangular views. Result is nonzero, where any of the inputs is
 * nonzero.
 * @param a first view
 * @param b second view
 * @return combined view
 */
inline TriangularViewCL combine(TriangularViewCL a, TriangularViewCL b) {
  typedef typename std::underlying_type<TriangularViewCL>::type underlying;
  return static_cast<TriangularViewCL>(static_cast<underlying>(a)
                                       | static_cast<underlying>(b));
}

/**
 * Determines common nonzero part of the inputs. Result is nonzero, where both
 * inputs are nonzero.
 * @param a first view
 * @param b second view
 * @return common nonzero part
 */
inline TriangularViewCL commonNonzeroPart(TriangularViewCL a,
                                          TriangularViewCL b) {
  typedef typename std::underlying_type<TriangularViewCL>::type underlying;
  return static_cast<TriangularViewCL>(static_cast<underlying>(a)
                                       & static_cast<underlying>(b));
}

/**
 * Check whether a view contains certain nonzero part
 * @param a view to check
 * @param b part to check for (usually `Lower` or `Upper`)
 * @return true, if `a` has part `b` nonzero
 */
inline bool containsNonzeroPart(TriangularViewCL a, TriangularViewCL b) {
  return static_cast<bool>(commonNonzeroPart(a, b));
}

/**
 * Transposes a triangular view - swaps lower and upper parts.
 * @param a view to transpose
 * @return transposition of input
 */
inline TriangularViewCL transpose(TriangularViewCL a) {
  if(a==TriangularViewCL::Lower){
    return TriangularViewCL::Upper;
  }
  if(a==TriangularViewCL::Upper){
    return TriangularViewCL::Lower;
  }
  return a;
}

/**
 * Inverts a triangular view. Parts that are zero in the input become nonzero in
 * output and vice versa.
 * @param a view to invert
 * @return inverted view
 */
inline TriangularViewCL invert(TriangularViewCL a) {
  typedef typename std::underlying_type<TriangularViewCL>::type underlying;
  return static_cast<TriangularViewCL>(static_cast<underlying>(TriangularViewCL::Entire) & ~static_cast<underlying>(a));
}

/**
 * Creates a triangular view from `Eigen::UpLoType`. `Eigen::Lower`,
 * `Eigen::StrictlyLower` and `Eigen::UnitLower` become
 * `TriangularViewCL::Lower`. Similar for `Upper`. Any other view becomes
 * `TriangularViewCL::Entire`.
 * @param a `UpLoType` to create a view from
 * @return triangular view
 */
inline TriangularViewCL fromEigenUpLoType(Eigen::UpLoType a) {
  if (a & Eigen::Lower) {
    return TriangularViewCL::Lower;
  }
  if (a & Eigen::Upper) {
    return TriangularViewCL::Upper;
  }
  return TriangularViewCL::Entire;
}

enum class TriangularMapCL { UpperToLower = 0, LowerToUpper = 1 };

// \cond
static const char *triangular_kernel_helpers = STRINGIFY(
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
