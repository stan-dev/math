#ifndef STAN_MATH_REV_CORE_PROFILING_HPP
#define STAN_MATH_REV_CORE_PROFILING_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <iostream>

namespace stan {
namespace math {

namespace internal {

inline auto profile_start(std::string name, profile_map& profiles) {
  profile_map::iterator it = profiles.find(name);
  if (it == profiles.end()) {
    profiles[name] = {};
  }
  if (profiles[name].meta.fwd_pass_active
      || profiles[name].meta.rev_pass_active) {
    throw std::runtime_error("FAIL");
  }
  profiles[name].meta.fwd_pass_start = std::chrono::steady_clock::now();
  profiles[name].meta.fwd_pass_active = true;
  reverse_pass_callback([name, &profiles]() mutable {
    profiles[name].rev_pass
        += std::chrono::duration<double>(std::chrono::steady_clock::now()
                                         - profiles[name].meta.rev_pass_start)
               .count();
    profiles[name].meta.rev_pass_active = false;
  });
}

inline auto profile_stop(std::string name, profile_map& profiles) {
  profile_map::iterator it = profiles.find(name);
  if (it == profiles.end()) {
    throw std::runtime_error("FAIL");
  }
  if (!profiles[name].meta.fwd_pass_active) {
    throw std::runtime_error("FAIL");
  }
  if (profiles[name].meta.rev_pass_active) {
    throw std::runtime_error("FAIL");
  }
  profiles[name].fwd_pass
      += std::chrono::duration<double>(std::chrono::steady_clock::now()
                                       - profiles[name].meta.fwd_pass_start)
             .count();
  profiles[name].meta.fwd_pass_active = false;

  reverse_pass_callback([name, &profiles]() mutable {
    profiles[name].meta.rev_pass_start = std::chrono::steady_clock::now();
    profiles[name].meta.rev_pass_active = true;
    profiles[name].count_rev++;
  });
}
}  // namespace internal

class profile {
  std::string name_;
  profile_map* profiles_;

 public:
  profile(std::string name, profile_map& profiles)
      : name_(name), profiles_(&profiles) {
    internal::profile_start(name_, *profiles_);
  }
  ~profile() { internal::profile_stop(name_, *profiles_); }
};

}  // namespace math
}  // namespace stan
#endif
