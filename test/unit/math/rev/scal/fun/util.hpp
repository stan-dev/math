#ifndef TEST_UNIT_MATH_REV_SCAL_FUN_UTIL_HPP
#define TEST_UNIT_MATH_REV_SCAL_FUN_UTIL_HPP

#include <stan/math/rev/scal.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <vector>

typedef stan::math::var AVAR;
typedef std::vector<AVAR> AVEC;
typedef std::vector<double> VEC;

AVEC createAVEC(AVAR x) {
  AVEC v;
  v.push_back(x);
  return v;
}
AVEC createAVEC(AVAR x1, AVAR x2) {
  AVEC v;
  v.push_back(x1);
  v.push_back(x2);
  return v;
}
AVEC createAVEC(AVAR x1, AVAR x2, AVAR x3) {
  AVEC v;
  v.push_back(x1);
  v.push_back(x2);
  v.push_back(x3);
  return v;
}
AVEC createAVEC(AVAR x1, AVAR x2, AVAR x3, AVAR x4) {
  AVEC v = createAVEC(x1, x2, x3);
  v.push_back(x4);
  return v;
}
AVEC createAVEC(AVAR x1, AVAR x2, AVAR x3, AVAR x4, AVAR x5) {
  AVEC v = createAVEC(x1, x2, x3, x4);
  v.push_back(x5);
  return v;
}
AVEC createAVEC(AVAR x1, AVAR x2, AVAR x3, AVAR x4, AVAR x5, AVAR x6) {
  AVEC v = createAVEC(x1, x2, x3, x4, x5);
  v.push_back(x6);
  return v;
}
AVEC createAVEC(AVAR x1, AVAR x2, AVAR x3, AVAR x4, AVAR x5, AVAR x6, AVAR x7) {
  AVEC v = createAVEC(x1, x2, x3, x4, x5, x6);
  v.push_back(x7);
  return v;
}
AVEC createAVEC(AVAR x1, AVAR x2, AVAR x3, AVAR x4, AVAR x5, AVAR x6, AVAR x7,
                AVAR x8) {
  AVEC v = createAVEC(x1, x2, x3, x4, x5, x6, x7);
  v.push_back(x8);
  return v;
}

VEC cgrad(AVAR f, AVAR x1) {
  AVEC x = createAVEC(x1);
  VEC g;
  f.grad(x, g);
  return g;
}

VEC cgrad(AVAR f, AVAR x1, AVAR x2) {
  AVEC x = createAVEC(x1, x2);
  VEC g;
  f.grad(x, g);
  return g;
}

VEC cgrad(AVAR f, AVAR x1, AVAR x2, AVAR x3) {
  AVEC x = createAVEC(x1, x2, x3);
  VEC g;
  f.grad(x, g);
  return g;
}

VEC cgrad(AVAR f, AVAR x1, AVAR x2, AVAR x3, AVAR x4) {
  AVEC x = createAVEC(x1, x2, x3, x4);
  VEC g;
  f.grad(x, g);
  return g;
}
#endif
