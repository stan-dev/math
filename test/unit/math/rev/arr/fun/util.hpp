#ifndef TEST_UNIT_MATH_REV_ARR_FUN_UTIL_HPP
#define TEST_UNIT_MATH_REV_ARR_FUN_UTIL_HPP

#include <vector>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/rev/arr/util.hpp>

typedef std::vector<AVAR> AVEC;
typedef std::vector<double> VEC;

VEC cgradvec(AVAR f, AVEC x) {
  VEC g;
  f.grad(x,g);
  return g;
}

#endif
