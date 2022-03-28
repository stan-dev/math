#ifndef STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/core.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <tuple>
#include <memory>
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
struct reduce_sum_impl<ReduceFunction, require_var_t<ReturnType>, ReturnType,
                       Vec, Args...> {
  struct scoped_args_tuple {
    ScopedChainableStack stack_;
    using args_tuple_t
        = std::tuple<decltype(deep_copy_vars(std::declval<Args>()))...>;
    std::unique_ptr<args_tuple_t> args_tuple_holder_;

    scoped_args_tuple() : stack_(), args_tuple_holder_(nullptr) {}
  };

  /**
   * This struct is used by the TBB to accumulate partial
   *  sums over consecutive ranges of the input. To distribute the workload,
   *  the TBB can split larger partial sums into smaller ones in which
   *  case the splitting copy constructor is used. It is designed to
   *  meet the Imperative form requirements of `tbb::parallel_reduce`.
   *
   * @note see link [here](https://tinyurl.com/vp7xw2t) for requirements.
   */
  struct recursive_reducer {
    const size_t num_vars_per_term_;
    const size_t num_vars_shared_terms_;  // Number of vars in shared arguments
    double* sliced_partials_;  // Points to adjoints of the partial calculations
    Vec vmapped_;
    std::stringstream msgs_;
    std::tuple<Args...> args_tuple_;
    scoped_args_tuple local_args_tuple_scope_;
    double sum_{0.0};
    Eigen::VectorXd args_adjoints_{0};

    template <typename VecT, typename... ArgsT>
    recursive_reducer(size_t num_vars_per_term, size_t num_vars_shared_terms,
                      double* sliced_partials, VecT&& vmapped, ArgsT&&... args)
        : num_vars_per_term_(num_vars_per_term),
          num_vars_shared_terms_(num_vars_shared_terms),
          sliced_partials_(sliced_partials),
          vmapped_(std::forward<VecT>(vmapped)),
          local_args_tuple_scope_(),
          args_tuple_(std::forward<ArgsT>(args)...) {}

    /*
     * This is the copy operator as required for tbb::parallel_reduce
     *   Imperative form. This requires sum_ and arg_adjoints_ be reset
     *   to zero since the newly created reducer is used to accumulate
     *   an independent partial sum.
     */
    recursive_reducer(recursive_reducer& other, tbb::split)
        : num_vars_per_term_(other.num_vars_per_term_),
          num_vars_shared_terms_(other.num_vars_shared_terms_),
          sliced_partials_(other.sliced_partials_),
          vmapped_(other.vmapped_),
          local_args_tuple_scope_(),
          args_tuple_(other.args_tuple_) {}

    /**
     * Compute, using nested autodiff, the value and Jacobian of
     *  `ReduceFunction` called over the range defined by r and accumulate those
     *  in member variable sum_ (for the value) and args_adjoints_ (for the
     *  Jacobian). The nested autodiff uses deep copies of the involved operands
     *  ensuring that no side effects are implied to the adjoints of the input
     *  operands which reside potentially on a autodiff tape stored in a
     *  different thread other than the current thread of execution. This
     * function may be called multiple times per object instantiation (so the
     * sum_ and args_adjoints_ must be accumulated, not just assigned).
     *
     * @param r Range over which to compute reduce_sum
     */
    inline void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty()) {
        return;
      }

      if (args_adjoints_.size() == 0) {
        args_adjoints_ = Eigen::VectorXd::Zero(num_vars_shared_terms_);
      }

      // Obtain reference to a local copy of all shared arguments that do
      // not point
      //   back to main autodiff stack

      if (!local_args_tuple_scope_.args_tuple_holder_) {
        // shared arguments need to be copied to reducer-specific
        // scope. In this case no need for zeroing adjoints, since the
        // fresh copy has all adjoints set to zero.
        local_args_tuple_scope_.stack_.execute([&]() {
          math::apply(
              [&](auto&&... args) {
                local_args_tuple_scope_.args_tuple_holder_ = std::make_unique<
                    typename scoped_args_tuple::args_tuple_t>(
                    deep_copy_vars(args)...);
              },
              args_tuple_);
        });
      } else {
        // set adjoints of shared arguments to zero
        local_args_tuple_scope_.stack_.execute([] { set_zero_all_adjoints(); });
      }

      auto& args_tuple_local = *(local_args_tuple_scope_.args_tuple_holder_);

      // Initialize nested autodiff stack
      const nested_rev_autodiff begin_nest;

      // Create nested autodiff copies of sliced argument that do not point
      //   back to main autodiff stack
      std::decay_t<Vec> local_sub_slice;
      local_sub_slice.reserve(r.size());
      for (size_t i = r.begin(); i < r.end(); ++i) {
        local_sub_slice.emplace_back(deep_copy_vars(vmapped_[i]));
      }

      // Perform calculation
      var sub_sum_v = math::apply(
          [&](auto&&... args) {
            return ReduceFunction()(local_sub_slice, r.begin(), r.end() - 1,
                                    &msgs_, args...);
          },
          args_tuple_local);

      // Compute Jacobian
      sub_sum_v.grad();

      // Accumulate value of reduce_sum
      sum_ += sub_sum_v.val();

      // Accumulate adjoints of sliced_arguments
      accumulate_adjoints(sliced_partials_ + r.begin() * num_vars_per_term_,
                          std::move(local_sub_slice));

      // Accumulate adjoints of shared_arguments
      math::apply(
          [&](auto&&... args) {
            accumulate_adjoints(args_adjoints_.data(), args...);
          },
          args_tuple_local);
    }

    /**
     * Join reducers. Accumuluate the value (sum_) and Jacobian (arg_adoints_)
     *   of the other reducer.
     *
     * @param rhs Another partial sum
     */
    inline void join(const recursive_reducer& rhs) {
      sum_ += rhs.sum_;
      if (args_adjoints_.size() != 0 && rhs.args_adjoints_.size() != 0) {
        args_adjoints_ += rhs.args_adjoints_;
      } else if (args_adjoints_.size() == 0 && rhs.args_adjoints_.size() != 0) {
        args_adjoints_ = rhs.args_adjoints_;
      }
      msgs_ << rhs.msgs_.str();
    }
  };

  /**
   * Call an instance of the function `ReduceFunction` on every element
   *   of an input sequence and sum these terms.
   *
   * This specialization is parallelized using tbb and works for reverse
   *   mode autodiff.
   *
   * ReduceFunction must define an operator() with the same signature as:
   *   var f(Vec&& vmapped_subset, int start, int end, std::ostream* msgs,
   * Args&&... args)
   *
   * `ReduceFunction` must be default constructible without any arguments
   *
   * Each call to `ReduceFunction` is responsible for computing the
   *   start through end (inclusive) terms of the overall sum. All args are
   * passed from this function through to the `ReduceFunction` instances.
   *   However, only the start through end (inclusive) elements of the vmapped
   * argument are passed to the `ReduceFunction` instances (as the
   * `vmapped_subset` argument).
   *
   * This function distributes computation of the desired sum and the Jacobian
   * of that sum over multiple threads by coordinating calls to `ReduceFunction`
   * instances. Results are stored as precomputed varis in the autodiff tree.
   *
   * If auto partitioning is true, break work into pieces automatically,
   *  taking grainsize as a recommended work size. The partitioning is
   *  not deterministic nor is the order guaranteed in which partial
   *  sums are accumulated. Due to floating point imprecisions this will likely
   *  lead to slight differences in the accumulated results between
   *  multiple runs. If false, break work deterministically into pieces smaller
   *  than or equal to grainsize and accumulate all the partial sums
   *  in the same order. This still may not achieve bitwise reproducibility.
   *
   * @param vmapped Vector containing one element per term of sum
   * @param auto_partitioning Work partitioning style
   * @param grainsize Suggested grainsize for tbb
   * @param[in, out] msgs The print stream for warning messages
   * @param args Shared arguments used in every sum term
   * @return Summation of all terms
   */
  inline var operator()(Vec&& vmapped, bool auto_partitioning, int grainsize,
                        std::ostream* msgs, Args&&... args) const {
    if (vmapped.empty()) {
      return var(0.0);
    }

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

    for (size_t i = 0; i < num_vars_sliced_terms; ++i) {
      partials[i] = 0.0;
    }

    recursive_reducer worker(num_vars_per_term, num_vars_shared_terms, partials,
                             std::forward<Vec>(vmapped),
                             std::forward<Args>(args)...);

    // we must use task isolation as described here:
    // https://software.intel.com/content/www/us/en/develop/documentation/tbb-documentation/top/intel-threading-building-blocks-developer-guide/task-isolation.html
    // this is to ensure that the thread local AD tape ressource is
    // not being modified from a different task which may happen
    // whenever this function is being used itself in a parallel
    // context (like running multiple chains for Stan)
    tbb::this_task_arena::isolate([&] {
      if (auto_partitioning) {
        tbb::parallel_reduce(
            tbb::blocked_range<std::size_t>(0, num_terms, grainsize), worker);
      } else {
        tbb::simple_partitioner partitioner;
        tbb::parallel_deterministic_reduce(
            tbb::blocked_range<std::size_t>(0, num_terms, grainsize), worker,
            partitioner);
      }
    });

    for (size_t i = 0; i < num_vars_shared_terms; ++i) {
      partials[num_vars_sliced_terms + i] = worker.args_adjoints_.coeff(i);
    }

    if (msgs) {
      *msgs << worker.msgs_.str();
    }

    return var(new precomputed_gradients_vari(
        worker.sum_, num_vars_sliced_terms + num_vars_shared_terms, varis,
        partials));
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
