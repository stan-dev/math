#ifndef STAN_MATH_PRIM_SCAL_META_LENGTH_MVT_HPP
#define STAN_MATH_PRIM_SCAL_META_LENGTH_MVT_HPP

#include <stdexcept>

namespace stan {

  template <typename T>
  size_t length_mvt(const T& ) {
    throw std::out_of_range("length_mvt passed to an unrecognized type.");
    return 1U;
  }

}
#endif

