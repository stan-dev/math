#ifndef STAN_MATH_REV_CORE_PROFILING_HPP
#define STAN_MATH_REV_CORE_PROFILING_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/core/profiling.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <iostream>
#include <sstream>
#include <thread>

namespace stan {
namespace math {

namespace internal {
template <typename T, require_all_var_t<T>* = nullptr>
inline auto profile_start(std::string name, profile_map& profiles) {
  profile_key key = {name, std::this_thread::get_id()};
  profile_map::iterator p = profiles.find(key);
  if (p == profiles.end()) {
    std::pair<profile_map::iterator, bool> new_profile
        = profiles.emplace(std::make_pair(key, profile_info({})));
    p = new_profile.first;
  }
  if (p->second.meta.fwd_pass_active) {
    std::ostringstream msg;
    msg << "Profile '" << name << "' ";
    msg << "already started!";
    throw std::runtime_error(std::string(msg.str()));
  }
  if (p->second.meta.rev_pass_active) {
    std::ostringstream msg;
    msg << "Forward and reverse pass active for profile '" << name << "'.";
    msg << "Please file a bug report!";
    throw std::runtime_error(std::string(msg.str()));
  }
  p->second.meta.fwd_pass_start = std::chrono::steady_clock::now();
  p->second.meta.fwd_pass_active = true;
  reverse_pass_callback([name, p]() mutable {
    profile_key key = {name, std::this_thread::get_id()};
    p->second.rev_pass
        += std::chrono::duration<double>(std::chrono::steady_clock::now()
                                         - p->second.meta.rev_pass_start)
               .count();
    p->second.meta.rev_pass_active = false;
  });
}

template <typename T, require_all_var_t<T>* = nullptr>
inline auto profile_stop(std::string name, profile_map& profiles) {
  profile_key key = {name, std::this_thread::get_id()};
  profile_map::iterator p = profiles.find(key);
  if (p == profiles.end()) {
    std::ostringstream msg;
    msg << "Stopping a non-registered profile '" << name
        << "'. Please file a bug report!";
    throw std::runtime_error("Stopping ");
  }
  if (!p->second.meta.fwd_pass_active) {
    std::ostringstream msg;
    msg << "Stopping forward pass profile '" << name
        << "' that was not started.";
    throw std::runtime_error(std::string(msg.str()));
  }
  p->second.fwd_pass
      += std::chrono::duration<double>(std::chrono::steady_clock::now()
                                       - p->second.meta.fwd_pass_start)
             .count();
  p->second.meta.fwd_pass_active = false;

  reverse_pass_callback([name, p]() mutable {
    p->second.meta.rev_pass_start = std::chrono::steady_clock::now();
    p->second.meta.rev_pass_active = true;
  });
}
}  // namespace internal

template <typename T>
class profile {
  std::string name_;
  profile_map* profiles_;

 public:
  profile(std::string name, profile_map& profiles)
      : name_(name), profiles_(&profiles) {
    internal::profile_start<T>(name_, *profiles_);
  }
  ~profile() { internal::profile_stop<T>(name_, *profiles_); }
};

}  // namespace math
}  // namespace stan
#endif
