#ifndef STAN_MATH_REV_CORE_PROFILING_HPP
#define STAN_MATH_REV_CORE_PROFILING_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <tbb/concurrent_unordered_map.h>
#include <iostream>
#include <sstream>
#include <thread>

namespace stan {
namespace math {

/**
 * Class used for storing profiling information.
 */
class profile_info {
 private:
  bool active_;

  double fwd_pass_time_;
  double rev_pass_time_;
  size_t n_fwd_AD_passes_;
  size_t n_fwd_no_AD_passes_;
  size_t n_rev_passes_;
  size_t chain_stack_size_sum_;
  size_t nochain_stack_size_sum_;
  std::chrono::time_point<std::chrono::steady_clock> fwd_pass_tp_;
  std::chrono::time_point<std::chrono::steady_clock> rev_pass_tp_;
  size_t start_chain_stack_size_;
  size_t start_nochain_stack_size_;

 public:
  profile_info()
      : active_(false),
        fwd_pass_time_(0.0),
        rev_pass_time_(0.0),
        n_fwd_AD_passes_(0),
        n_fwd_no_AD_passes_(0),
        n_rev_passes_(0),
        chain_stack_size_sum_(0),
        nochain_stack_size_sum_(0),
        fwd_pass_tp_(std::chrono::steady_clock::now()),
        rev_pass_tp_(std::chrono::steady_clock::now()),
        start_chain_stack_size_(0),
        start_nochain_stack_size_(0) {}

  bool is_active() const noexcept { return active_; }

  template <typename T>
  void fwd_pass_start() {
    if (!is_constant<T>::value) {
      start_chain_stack_size_ = ChainableStack::instance_->var_stack_.size();
      start_nochain_stack_size_
          = ChainableStack::instance_->var_nochain_stack_.size();
    }
    fwd_pass_tp_ = std::chrono::steady_clock::now();
    active_ = true;
  }

  template <typename T>
  void fwd_pass_stop() {
    if (!is_constant<T>::value) {
      n_fwd_AD_passes_++;
      chain_stack_size_sum_ += (ChainableStack::instance_->var_stack_.size()
                                - start_chain_stack_size_ - 1);
      nochain_stack_size_sum_
          += (ChainableStack::instance_->var_nochain_stack_.size()
              - start_nochain_stack_size_);
    } else {
      n_fwd_no_AD_passes_++;
    }
    fwd_pass_time_ += std::chrono::duration<double>(
                          std::chrono::steady_clock::now() - fwd_pass_tp_)
                          .count();
    active_ = false;
  }

  void rev_pass_start() { rev_pass_tp_ = std::chrono::steady_clock::now(); }

  void rev_pass_stop() {
    rev_pass_time_ += std::chrono::duration<double>(
                          std::chrono::steady_clock::now() - rev_pass_tp_)
                          .count();
    n_rev_passes_++;
  }

  size_t get_chain_stack_used() const noexcept {
    return chain_stack_size_sum_;
  };

  size_t get_nochain_stack_used() const noexcept {
    return nochain_stack_size_sum_;
  };

  size_t get_num_no_AD_fwd_passes() const noexcept {
    return n_fwd_no_AD_passes_;
  }

  size_t get_num_AD_fwd_passes() const noexcept { return n_fwd_AD_passes_; }

  size_t get_num_fwd_passes() const noexcept {
    return n_fwd_AD_passes_ + n_fwd_no_AD_passes_;
  }

  double get_fwd_time() const noexcept { return fwd_pass_time_; }

  size_t get_num_rev_passes() const noexcept { return n_rev_passes_; }

  double get_rev_time() const noexcept { return rev_pass_time_; }
};

using profile_key = std::pair<std::string, std::thread::id>;

namespace internal {
struct hash_profile_key {
  std::size_t operator()(const profile_key& key) const {
    return std::hash<std::string>()(key.first)
           ^ std::hash<std::thread::id>()(key.second);
  }
};
struct equal_profile_key {
  bool operator()(const profile_key& lhs, const profile_key& rhs) const {
    return lhs.first == rhs.first && lhs.second == rhs.second;
  }
};

}  // namespace internal

using profile_map = tbb::concurrent_unordered_map<profile_key, profile_info,
                                                  internal::hash_profile_key,
                                                  internal::equal_profile_key>;

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
  profile_key key_;
  profile_info* profile_;

 public:
  profile(std::string name, profile_map& profiles)
      : key_({name, std::this_thread::get_id()}) {
    profile_map::iterator p = profiles.find(key_);
    if (p == profiles.end()) {
      profiles[key_] = profile_info();
    }
    profile_ = &profiles[key_];
    if (profile_->is_active()) {
      std::ostringstream msg;
      msg << "Profile '" << key_.first << "' already started!";
      throw std::runtime_error(msg.str());
    }
    profile_->fwd_pass_start<T>();
    if (!is_constant<T>::value) {
      reverse_pass_callback(
          [profile = this->profile_]() mutable { profile->rev_pass_stop(); });
    }
  }
  ~profile() {
    profile_->fwd_pass_stop<T>();
    if (!is_constant<T>::value) {
      reverse_pass_callback(
          [profile = this->profile_]() mutable { profile->rev_pass_start(); });
    }
  }
};

}  // namespace math
}  // namespace stan
#endif
