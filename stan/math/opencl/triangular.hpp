#ifndef STAN_MATH_OPENCL_TRIANGULAR_HPP
#define STAN_MATH_OPENCL_TRIANGULAR_HPP
#ifdef STAN_OPENCL
namespace stan {
namespace math {
enum class TriangularViewCL { Diagonal = 0, Lower = 1, Upper = 2, Entire = 3 };

TriangularViewCL operator| (TriangularViewCL a, TriangularViewCL b){
  typedef typename std::underlying_type<TriangularViewCL>::type underlying;
  return static_cast<TriangularViewCL>(static_cast<underlying>(a) | static_cast<underlying>(b));
}

TriangularViewCL operator|= (TriangularViewCL& a, TriangularViewCL b){
  a = a | b;
  return a;
}

TriangularViewCL operator& (TriangularViewCL a, TriangularViewCL b){
  typedef typename std::underlying_type<TriangularViewCL>::type underlying;
  return static_cast<TriangularViewCL>(static_cast<underlying>(a) & static_cast<underlying>(b));
}

TriangularViewCL operator&= (TriangularViewCL& a, TriangularViewCL b){
  a = a & b;
  return a;
}

TriangularViewCL transpose(TriangularViewCL a){
  switch (a){
    case TriangularViewCL::Lower:
      return TriangularViewCL::Upper;
    case TriangularViewCL::Upper:
      return TriangularViewCL::Lower;
    default:
      return a;
  }
}

TriangularViewCL invert(TriangularViewCL a){
  switch (a){
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

enum class TriangularMapCL { UpperToLower = 0, LowerToUpper = 1 };
}  // namespace math
}  // namespace stan
#endif
#endif
