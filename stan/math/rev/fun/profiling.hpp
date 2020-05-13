#ifndef STAN_MATH_REV_FUN_PROFILING_HPP
#define STAN_MATH_REV_FUN_PROFILING_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/acosh.hpp>
#include <stan/math/prim/fun/isnan.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/arg.hpp>
#include <stan/math/rev/fun/cosh.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/polar.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

void eval() { } 

template <typename T, typename = require_not_eigen_t<T>>
void eval(T a) { } 

template <typename T, require_eigen_t<T>* = nullptr>
void eval(T a) { a.eval(); } 

template <typename T, typename... Types> 
void eval(T var1, Types... var2) 
{ 
    eval(var1);
  
    eval(var2...) ; 
}

namespace internal {
 template <typename... Types> 
class start_profiling_vari : public vari {
 int id_;
 profiles& pp;

 public:
  start_profiling_vari(int id, profiles& p, Types... args) : vari(0), pp(p) {
      eval(args...);
      id_ = id;
      std::cout << "Forward pass start: " << id << std::endl;

  }
  void chain() {
    std::cout << "Reverse pass end: " << id_ << std::endl;
  }
};

 template <typename... Types> 
class stop_profiling_vari : public vari {
 int id_;
 profiles& pp;
 public:
  stop_profiling_vari(int id, profiles& p, Types... args) : vari(0), pp(p) {
      eval(args...);
      id_ = id;
      std::cout << "Forward pass stop: " << id << std::endl;

  }
  void chain() {
    //eval(args..);
    std::cout << "Reverse pass start: " << id_ << std::endl;
  }
};
}  // namespace internal

template <typename... Types> 
inline var start_profiling(int id, profiles& p, Types... args) {
  return var(new internal::start_profiling_vari(id, p, args...));
}

template <typename... Types> 
inline var stop_profiling(int id, profiles& p, Types... args) {
  return var(new internal::stop_profiling_vari(id, p, args...));
}

}  // namespace math
}  // namespace stan
#endif
