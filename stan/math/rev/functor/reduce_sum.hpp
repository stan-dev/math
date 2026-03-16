#ifndef STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/team_thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename ReduceFunction, typename ReturnType, typename Vec,
          typename... Args>
struct reduce_sum_impl<ReduceFunction, require_var_t<ReturnType>, ReturnType,
                       Vec, Args...> {
  struct scoped_args_tuple {
    ScopedChainableStack stack_;
    using args_tuple_t = std::tuple<decltype(
        deep_copy_vars(std::declval<std::decay_t<Args>>()))...>;
    std::unique_ptr<args_tuple_t> args_tuple_holder_;
    scoped_args_tuple() : stack_(), args_tuple_holder_(nullptr) {}
  };

  struct recursive_reducer {
    using VecRef = std::decay_t<Vec>;
    using args_ptrs_t = std::tuple<const std::decay_t<Args>*...>;

    // Apply a callable to tuple of pointers, dereferencing each pointer.
    template <typename Fn, typename Tuple, std::size_t... I>
    static inline decltype(auto) apply_ptr_tuple_impl(
        Fn&& fn, Tuple& t, std::index_sequence<I...>) {
      return std::forward<Fn>(fn)((*std::get<I>(t))...);
    }
    template <typename Fn, typename... Ts>
    static inline decltype(auto) apply_ptr_tuple(Fn&& fn,
                                                 std::tuple<const Ts*...>& t) {
      return apply_ptr_tuple_impl(std::forward<Fn>(fn), t,
                                  std::index_sequence_for<Ts...>{});
    }

    const std::size_t num_vars_per_term_;
    const std::size_t num_vars_shared_terms_;
    double* sliced_partials_;

    const VecRef* vmapped_;
    args_ptrs_t args_ptrs_;

    // msgs only if requested
    std::unique_ptr<std::stringstream> msgs_;
    std::ostream* msgs_out_;

    scoped_args_tuple local_args_tuple_scope_;

    double sum_{0.0};
    std::vector<double> args_adjoints_;  // faster than Eigen::VectorXd

    // Reusable buffer to avoid realloc per chunk
    VecRef local_sub_slice_;
    std::size_t reserved_ = 0;

    recursive_reducer(std::size_t num_vars_per_term,
                      std::size_t num_vars_shared_terms,
                      double* sliced_partials, const VecRef& vmapped,
                      std::ostream* msgs, const std::decay_t<Args>&... args)
        : num_vars_per_term_(num_vars_per_term),
          num_vars_shared_terms_(num_vars_shared_terms),
          sliced_partials_(sliced_partials),
          vmapped_(&vmapped),
          args_ptrs_(&args...),
          msgs_out_(msgs) {
      if (msgs_out_)
        msgs_ = std::make_unique<std::stringstream>();
    }

    inline void operator()(std::size_t begin, std::size_t end) {
      if (begin >= end)
        return;

      if (args_adjoints_.empty()) {
        args_adjoints_.assign(num_vars_shared_terms_, 0.0);
      }

      if (!local_args_tuple_scope_.args_tuple_holder_) {
        local_args_tuple_scope_.stack_.execute([&]() {
          apply_ptr_tuple(
              [&](auto const&... a) {
                local_args_tuple_scope_.args_tuple_holder_ = std::make_unique<
                    typename scoped_args_tuple::args_tuple_t>(
                    deep_copy_vars(a)...);
              },
              args_ptrs_);
        });
      } else {
        local_args_tuple_scope_.stack_.execute([] { set_zero_all_adjoints(); });
      }

      auto& args_tuple_local = *(local_args_tuple_scope_.args_tuple_holder_);

      const nested_rev_autodiff begin_nest;

      // Reuse per-worker buffer
      const std::size_t n = end - begin;
      local_sub_slice_.clear();
      if (reserved_ < n) {
        local_sub_slice_.reserve(n);
        reserved_ = n;
      }

      for (std::size_t i = begin; i < end; ++i) {
        local_sub_slice_.emplace_back(deep_copy_vars((*vmapped_)[i]));
      }

      std::ostream* local_msgs
          = msgs_ ? static_cast<std::ostream*>(msgs_.get()) : nullptr;

      var sub_sum_v = math::apply(
          [&](auto&&... args_local) {
            return ReduceFunction()(local_sub_slice_, begin, end - 1,
                                    local_msgs, args_local...);
          },
          args_tuple_local);

      sub_sum_v.grad();
      sum_ += sub_sum_v.val();

      accumulate_adjoints(sliced_partials_ + begin * num_vars_per_term_,
                          std::move(local_sub_slice_));

      // local_sub_slice_ got moved-from; restore it to a valid empty state
      local_sub_slice_.clear();

      math::apply(
          [&](auto&&... args_local) {
            accumulate_adjoints(args_adjoints_.data(), args_local...);
          },
          args_tuple_local);
    }
  };

  inline var operator()(Vec&& vmapped, bool /*auto_partitioning*/,
                        int grainsize, std::ostream* msgs,
                        Args&&... args) const {
    if (vmapped.empty())
      return var(0.0);

    const std::size_t num_terms = vmapped.size();
    const std::size_t num_vars_per_term = count_vars(vmapped[0]);
    const std::size_t num_vars_sliced_terms = num_terms * num_vars_per_term;
    const std::size_t num_vars_shared_terms = count_vars(args...);

    vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(
        num_vars_sliced_terms + num_vars_shared_terms);
    double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
        num_vars_sliced_terms + num_vars_shared_terms);

    save_varis(varis, vmapped);
    save_varis(varis + num_vars_sliced_terms, args...);

    for (std::size_t i = 0; i < num_vars_sliced_terms; ++i)
      partials[i] = 0.0;

    auto& pool = stan::math::TeamThreadPool::instance();
    const std::size_t max_team = pool.team_size();

    // Choose workers. (Caller participates, so total participants = n)
    std::size_t n
        = std::min<std::size_t>(max_team, num_terms == 0 ? 1 : num_terms);
    if (n < 1)
      n = 1;

    // Chunking: default to ~2 chunks per participant (lower overhead).
    std::size_t gs;
    if (grainsize > 0) {
      gs = static_cast<std::size_t>(grainsize);
      if (gs < 1)
        gs = 1;
    } else {
      const std::size_t target_chunks = n * 2;
      gs = (num_terms + target_chunks - 1) / target_chunks;
      if (gs < 1)
        gs = 1;
    }

    // Serial cutoff: if too few terms, don't parallelize.
    if (n == 1 || num_terms <= gs) {
      recursive_reducer r(num_vars_per_term, num_vars_shared_terms, partials,
                          vmapped, msgs, args...);
      r(0, num_terms);

      // write shared adjoints
      if (!r.args_adjoints_.empty()) {
        for (std::size_t i = 0; i < num_vars_shared_terms; ++i) {
          partials[num_vars_sliced_terms + i] = r.args_adjoints_[i];
        }
      } else {
        for (std::size_t i = 0; i < num_vars_shared_terms; ++i) {
          partials[num_vars_sliced_terms + i] = 0.0;
        }
      }

      if (msgs && r.msgs_)
        *msgs << r.msgs_->str();

      return var(new precomputed_gradients_vari(
          r.sum_, num_vars_sliced_terms + num_vars_shared_terms, varis,
          partials));
    }

    // One reducer per participant (0..n-1) for static partitioning.
    // NOTE: we avoid copying vmapped/args by taking references/pointers inside
    // reducer.
    std::vector<std::unique_ptr<recursive_reducer>> workers;
    workers.reserve(n);
    for (std::size_t tid = 0; tid < n; ++tid) {
      workers.emplace_back(std::make_unique<recursive_reducer>(
          num_vars_per_term, num_vars_shared_terms, partials, vmapped, msgs,
          args...));
    }
    /*
    std::cout <<
    "--------------------------------------------------------------------------------"
    << std::endl
              << "worker count = " << pool.worker_count() << std::endl
              << "team size = " << pool.team_size() << std::endl
              << "gs = " << gs << std::endl
              << std::endl << std::endl;
    */

    // Static partition: each participant gets a contiguous block once
    pool.parallel_region(n, [&](std::size_t tid) {
      const std::size_t b0 = (num_terms * tid) / n;
      const std::size_t b1 = (num_terms * (tid + 1)) / n;
      if (b0 < b1) {
        (*workers[tid])(b0, b1);
      }
    });

    // Aggregate
    double total_sum = 0.0;
    std::vector<double> shared_adj(num_vars_shared_terms, 0.0);
    std::stringstream all_msgs;

    for (auto& wptr : workers) {
      auto& w = *wptr;
      total_sum += w.sum_;
      if (!w.args_adjoints_.empty()) {
        for (std::size_t i = 0; i < num_vars_shared_terms; ++i) {
          shared_adj[i] += w.args_adjoints_[i];
        }
      }
      if (msgs && w.msgs_) {
        all_msgs << w.msgs_->str();
      }
    }

    for (std::size_t i = 0; i < num_vars_shared_terms; ++i) {
      partials[num_vars_sliced_terms + i] = shared_adj[i];
    }

    if (msgs)
      *msgs << all_msgs.str();

    return var(new precomputed_gradients_vari(
        total_sum, num_vars_sliced_terms + num_vars_shared_terms, varis,
        partials));
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
