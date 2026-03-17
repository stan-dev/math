#ifndef TEST_UNIT_MULTIPLE_TRANSLATION_UNITS_INCLUDES_HPP
#define TEST_UNIT_MULTIPLE_TRANSLATION_UNITS_INCLUDES_HPP

// Keep this list aligned with the public Stan Math umbrella headers. Each
// translation unit test source includes this header so linking catches ODR
// violations from any externally linked header-defined function reachable from
// the public entry points.
#include <stan/math.hpp>
#include <stan/math/fwd.hpp>
#include <stan/math/mix.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>

#endif
