#ifndef STAN_MATH_REV_CORE_PROFILING_HPP
#define STAN_MATH_REV_CORE_PROFILING_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
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
/**
 * Starts profiling the forward pass for the profile
 * with the specified name and adds a var that stops
 * profiling the reverse pass in the reverse callback
 * if the specified template parameter T is a var.
 *
 * @tparam T profile type which must not a var
 *
 * @param name the name of the profile to start
 * @param profile reference to the object storing profiling data
 */
template <typename T>
inline auto profile_start(std::string name, profile_info& profile) {
  using std::chrono::duration;
  using std::chrono::steady_clock;
  if (profile.meta.fwd_pass_active) {
    std::ostringstream msg;
    msg << "Profile '" << name << "' ";
    msg << "already started!";
    throw std::runtime_error(msg.str());
  }
  profile.meta.fwd_pass_start_tp = steady_clock::now();
  profile.meta.fwd_pass_active = true;
  if (!is_constant<T>::value) {
     profile.meta.start_chain_stack_size
      = ChainableStack::instance_->var_stack_.size();
    reverse_pass_callback([&profile]() mutable {
      profile.rev_pass_time += duration<double>(steady_clock::now()
                                                - profile.meta.rev_pass_start_tp)
                                  .count();
      profile.meta.rev_pass_active = false;
      profile.n_rev_pass++;
    });
  } 
}

/**
 * Stops profiling the forward pass for the profile
 * with the specified name and adds a var that starts
 * profiling the reverse pass in the reverse callback
 * if the specified template parameter T is a var.
 *
 * @tparam T profile type which must not a var
 *
 * @param name the name of the profile to stop
 * @param profile reference to the object storing profiling data
 */
template <typename T>
inline auto profile_stop(std::string name, profile_info& profile) {
  using std::chrono::duration;
  using std::chrono::steady_clock;
  profile.fwd_pass_time
      += duration<double>(steady_clock::now() - profile.meta.fwd_pass_start_tp)
             .count();
  profile.meta.fwd_pass_active = false;
  profile.n_fwd_pass++;
  if (!is_constant<T>::value) {
    profile.chain_stack_size_sum += (ChainableStack::instance_->var_stack_.size()
                                    - profile.meta.start_chain_stack_size - 1);
    profile.nochain_stack_size_sum
        += (ChainableStack::instance_->var_nochain_stack_.size()
            - profile.meta.start_nochain_stack_size);
    reverse_pass_callback([&profile]() mutable {
      profile.meta.rev_pass_start_tp = steady_clock::now();
      profile.meta.rev_pass_active = true;
    });
  }
}
}  // namespace internal

/**
 * Profiles C++ lines where the object is in scope.
 * When T is var, the constructor starts the profile for the forward pass
 * and places a var with a callback to stop the reverse pass on the
 * AD tape. The destructor stops the profile for the forward pass and
 * places a var with a callback to start the profile for the reverse pass.
 * When T is not var, the constructor and destructor only profile the
 *
 *
 * @tparam T type of profile class. If var, the created object is used
 * to profile reverse mode AD. Only profiles the forward pass otherwise.
 */
template <typename T>
class profile {
  std::string name_;
  profile_info* profile_;

 public:
  profile(std::string name, profile_map& profiles) : name_(name) {
    profile_key key = {name_, std::this_thread::get_id()};
    profile_map::iterator p = profiles.find(key);
    if (p == profiles.end()) {
      profiles[key] = {};
    }
    profile_ = &profiles[key];
    internal::profile_start<T>(name_, *profile_);
  }
  ~profile() { internal::profile_stop<T>(name_, *profile_); }
};

}  // namespace math
}  // namespace stan
#endif
