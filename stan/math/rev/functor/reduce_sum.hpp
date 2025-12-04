#ifndef STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/simple_thread_pool.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <memory>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

namespace stan {
namespace math {
namespace internal {

/**
 * Var specialization of reduce_sum_impl
 *
 * @tparam ReduceFunction Type of reducer function
 * @tparam ReturnType Must be var
 * @tparam Vec Type of sliced argument
 * @tparam Args Types of shared arguments
 */
template <typename ReduceFunction, typename ReturnType, typename Vec,
          typename... Args>
struct reduce_sum_impl<ReduceFunction, require_var_t<ReturnType>, ReturnType, Vec,
                       Args...> {
  struct scoped_args_tuple {
    ScopedChainableStack stack_;
    using args_tuple_t
        = std::tuple<decltype(deep_copy_vars(std::declval<Args>()))...>;
    std::unique_ptr<args_tuple_t> args_tuple_holder_;

    scoped_args_tuple() : stack_(), args_tuple_holder_(nullptr) {}
  };

  struct recursive_reducer {
    const std::size_t num_vars_per_term_;
    const std::size_t num_vars_shared_terms_;
    double* sliced_partials_;

    Vec vmapped_;
    std::stringstream msgs_;
    std::tuple<Args...> args_tuple_;
    scoped_args_tuple local_args_tuple_scope_;

    double sum_{0.0};
    Eigen::VectorXd args_adjoints_{0};

    template <typename VecT, typename... ArgsT>
    recursive_reducer(std::size_t num_vars_per_term,
                      std::size_t num_vars_shared_terms,
                      double* sliced_partials,
                      VecT&& vmapped,
                      ArgsT&&... args)
        : num_vars_per_term_(num_vars_per_term),
          num_vars_shared_terms_(num_vars_shared_terms),
          sliced_partials_(sliced_partials),
          vmapped_(std::forward<VecT>(vmapped)),
          args_tuple_(std::forward<ArgsT>(args)...) {}

    inline void operator()(std::size_t begin, std::size_t end) {
      if (begin >= end) {
        return;
      }

      if (args_adjoints_.size() == 0) {
        args_adjoints_ = Eigen::VectorXd::Zero(num_vars_shared_terms_);
      }

      // Obtain reference to a local copy of all shared arguments that do not
      // point back to main autodiff stack.
      if (!local_args_tuple_scope_.args_tuple_holder_) {
        // First time: copy shared arguments into reducer-scoped space.
        local_args_tuple_scope_.stack_.execute([&]() {
          math::apply(
              [&](auto&&... args) {
                local_args_tuple_scope_.args_tuple_holder_ =
                    std::make_unique<typename scoped_args_tuple::args_tuple_t>(
                        deep_copy_vars(args)...);
              },
              args_tuple_);
        });
      } else {
        // Subsequent calls: reset adjoints of shared arguments.
        local_args_tuple_scope_.stack_.execute([] { set_zero_all_adjoints(); });
      }

      auto& args_tuple_local = *(local_args_tuple_scope_.args_tuple_holder_);

      // Start a nested autodiff scope for this chunk.
      const nested_rev_autodiff begin_nest;

      // Copy sliced terms into this nested stack.
      std::decay_t<Vec> local_sub_slice;
      local_sub_slice.reserve(end - begin);
      for (std::size_t i = begin; i < end; ++i) {
        local_sub_slice.emplace_back(deep_copy_vars(vmapped_[i]));
      }

      // Run user reducer over [begin, end-1]
      var sub_sum_v = math::apply(
          [&](auto&&... args) {
            return ReduceFunction()(local_sub_slice, begin, end - 1, &msgs_,
                                    args...);
          },
          args_tuple_local);

      // Compute partial derivatives within this nested stack.
      sub_sum_v.grad();

      // Accumulate value.
      sum_ += sub_sum_v.val();

      // Accumulate adjoints of sliced arguments into the shared partials array.
      accumulate_adjoints(sliced_partials_ + begin * num_vars_per_term_,
                          std::move(local_sub_slice));

      // Accumulate adjoints of shared arguments into this reducer's buffer.
      math::apply(
          [&](auto&&... args) {
            accumulate_adjoints(args_adjoints_.data(), args...);
          },
          args_tuple_local);
    }
  };

  inline var operator()(Vec&& vmapped, bool /*auto_partitioning*/, int grainsize,
                        std::ostream* msgs, Args&&... args) const {
    if (vmapped.empty()) {
      return var(0.0);
    }

    const std::size_t num_terms = vmapped.size();
    const std::size_t num_vars_per_term = count_vars(vmapped[0]);
    const std::size_t num_vars_sliced_terms = num_terms * num_vars_per_term;
    const std::size_t num_vars_shared_terms = count_vars(args...);

    vari** varis
        = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
            num_vars_sliced_terms + num_vars_shared_terms);
    double* partials
        = ChainableStack::instance_->memalloc_.alloc_array<double>(
            num_vars_sliced_terms + num_vars_shared_terms);

    save_varis(varis, vmapped);
    save_varis(varis + num_vars_sliced_terms, args...);

    for (std::size_t i = 0; i < num_vars_sliced_terms; ++i) {
      partials[i] = 0.0;
    }

    // Decide chunk size (grainsize) and workers.
    auto& pool = stan::math::SimpleThreadPool::instance();
    std::size_t pool_threads = pool.thread_count();
    if (pool_threads == 0) pool_threads = 1;

    std::size_t gs;
    if (grainsize > 0) {
      gs = static_cast<std::size_t>(grainsize);
    } else {
      // Heuristic: target ~8 chunks per worker.
      const std::size_t target_chunks = pool_threads * 8;
      gs = (num_terms + target_chunks - 1) / target_chunks;
      if (gs < 1) gs = 1;
    }

    // Serial cutoffs.
    if (pool_threads == 1 || num_terms <= gs) {
      recursive_reducer r(num_vars_per_term, num_vars_shared_terms,
                          partials, vmapped, args...);
      r(0, num_terms);

      for (std::size_t i = 0; i < num_vars_shared_terms; ++i) {
        partials[num_vars_sliced_terms + i]
            = (r.args_adjoints_.size() != 0) ? r.args_adjoints_.coeff(i) : 0.0;
      }
      if (msgs) {
        *msgs << r.msgs_.str();
      }

      return var(new precomputed_gradients_vari(
          r.sum_, num_vars_sliced_terms + num_vars_shared_terms, varis,
          partials));
    }

    const std::size_t num_workers = std::min(pool_threads, num_terms);

    // One reducer per worker (reused across chunks). This lets each worker reuse
    // its deep-copied shared arguments and local stack.
    std::vector<std::unique_ptr<recursive_reducer>> workers;
    workers.reserve(num_workers);
    for (std::size_t t = 0; t < num_workers; ++t) {
      workers.emplace_back(std::make_unique<recursive_reducer>(
          num_vars_per_term, num_vars_shared_terms, partials, vmapped, args...));
    }

    std::atomic<std::size_t> next{0};

    // Dynamic scheduling: each worker pulls the next chunk.
    pool.parallel_region(num_workers, [&](std::size_t tid) {
      auto& w = *workers[tid];
      while (true) {
        const std::size_t begin =
            next.fetch_add(gs, std::memory_order_relaxed);
        if (begin >= num_terms) break;
        const std::size_t end = std::min<std::size_t>(begin + gs, num_terms);
        w(begin, end);
      }
    });

    // Aggregate results.
    double total_sum = 0.0;
    Eigen::VectorXd shared_adjoints =
        Eigen::VectorXd::Zero(num_vars_shared_terms);
    std::stringstream all_msgs;

    for (auto& wptr : workers) {
      auto& w = *wptr;
      total_sum += w.sum_;
      if (w.args_adjoints_.size() != 0) {
        shared_adjoints += w.args_adjoints_;
      }
      all_msgs << w.msgs_.str();
    }

    for (std::size_t i = 0; i < num_vars_shared_terms; ++i) {
      partials[num_vars_sliced_terms + i] = shared_adjoints.coeff(i);
    }

    if (msgs) {
      *msgs << all_msgs.str();
    }

    return var(new precomputed_gradients_vari(
        total_sum, num_vars_sliced_terms + num_vars_shared_terms, varis,
        partials));
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
