#ifndef STAN_MATH_OPENCL_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_TRIANGULAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/stringify.hpp>
#include <Eigen/Core>
#include <type_traits>

namespace stan {
namespace math {
enum class TriangularViewCL { Diagonal = 0, Lower = 1, Upper = 2, Entire = 3 };

inline TriangularViewCL combine(TriangularViewCL a, TriangularViewCL b) {
  typedef typename std::underlying_type<TriangularViewCL>::type underlying;
  return static_cast<TriangularViewCL>(static_cast<underlying>(a)
                                       | static_cast<underlying>(b));
}

inline TriangularViewCL commonNonzeroPart(TriangularViewCL a,
                                          TriangularViewCL b) {
  typedef typename std::underlying_type<TriangularViewCL>::type underlying;
  return static_cast<TriangularViewCL>(static_cast<underlying>(a)
                                       & static_cast<underlying>(b));
}

inline bool containsNonzeroPart(TriangularViewCL a, TriangularViewCL b) {
  return static_cast<bool>(commonNonzeroPart(a, b));
}

inline TriangularViewCL transpose(TriangularViewCL a) {
  switch (a) {
    case TriangularViewCL::Lower:
      return TriangularViewCL::Upper;
    case TriangularViewCL::Upper:
      return TriangularViewCL::Lower;
    default:
      return a;
  }
}

inline TriangularViewCL invert(TriangularViewCL a) {
  switch (a) {
    case TriangularViewCL::Lower:
      return TriangularViewCL::Upper;
    case TriangularViewCL::Upper:
      return TriangularViewCL::Lower;
    case TriangularViewCL::Entire:
      return TriangularViewCL::Diagonal;
    default:
      return TriangularViewCL::Entire;
  }
}

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
    int combine(int a, int b) { return a | b; }

    int commonNonzeroPart(int a, int b) { return a & b; }

    bool containsNonzeroPart(int a, int b) { return commonNonzeroPart(a, b); }
    // \cond
);
// \endcond

}  // namespace math
}  // namespace stan
#endif
#endif
