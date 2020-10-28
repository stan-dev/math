#ifndef STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/rev/core.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <tuple>
#include <memory>
#include <utility>
#include <vector>

#include <tbb/enumerable_thread_specific.h>

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

  using args_tuple_t = std::tuple<decltype(deep_copy_vars(std::declval<Args>()))...>;

  /**
   * This is a data structure that we use for the forward computation only!
   *
   * Each thread will get one copy of this data structure (it is stored in
   *   a threadlocal thing)
   *
   * Pieces of this data structure are passed in to the revese mode vari, but
   *   the data structure itself is not used there
   *
   * It holds:
   *  1. A pointer to the Stack on which we do our calculations
   *  2. A pointer to a local copy of the input arguments (copied to our Stack)
   *  3. A list of the last var computed in every work chunk, this is the
   *     result of one partial sum. Since each thread may do many partial
   *     sums, this needs to be dynamically sized
   */
  struct stuff {
    ScopedChainableStack* stack_;
    args_tuple_t* args_tuple_;
    std::vector<var> tails_;

    stuff(ScopedChainableStack* stack, Args&&... args) :
      stack_(stack) {
      stack_->execute([&]() {
	args_tuple_ = new args_tuple_t(deep_copy_vars(args)...);
      });
    }
  };

  using local_stuff = tbb::enumerable_thread_specific<stuff>;

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
    local_stuff& tlstuff_;
    Vec vmapped_;
    std::ostream* msgs_;
    double sum_;

    template <typename VecT, typename... ArgsT>
    recursive_reducer(local_stuff& tlstuff, VecT&& vmapped,
                      std::ostream* msgs)
        : tlstuff_(tlstuff),
          vmapped_(std::forward<VecT>(vmapped)),
          msgs_(msgs),
	  sum_(0.0) {}

    /*
     * This is the copy operator as required for tbb::parallel_reduce
     *   Imperative form. This requires sum_ and arg_adjoints_ be reset
     *   to zero since the newly created reducer is used to accumulate
     *   an independent partial sum.
     */
    recursive_reducer(recursive_reducer& other, tbb::split)
      : tlstuff_(other.tlstuff_),
	vmapped_(other.vmapped_),
	msgs_(other.msgs_),
	sum_(0.0) {}

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

      /**
       * Obtain reference to thread local copy of all shared arguments
       * 
       * This object is a `stuff` object from above
       */
      auto& tl = tlstuff_.local();

      // Do the forward pass calculation
      tl.stack_->execute([&] {
	std::decay_t<Vec> local_sub_slice;
        local_sub_slice.reserve(r.size());
        for (size_t i = r.begin(); i < r.end(); ++i) {
          local_sub_slice.emplace_back(vmapped_[i]);
        }

        // Perform calculation
        var sub_sum_v = apply([&](auto&&... args) {
	    return ReduceFunction()(local_sub_slice, r.begin(), r.end() - 1,
				    msgs_, args...);
	  }, *tl.args_tuple_);

	// Save the `var` from the output of each partial_sum
	tl.tails_.push_back(sub_sum_v);
	sum_ += sub_sum_v.val();
      });
    }

    /**
     * Join reducers. Accumuluate the value (sum_) of the other reducer.
     *
     * @param rhs Another partial sum
     */
    inline void join(const recursive_reducer& rhs) {
      sum_ += rhs.sum_;
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
    const std::size_t num_terms = vmapped.size();

    if (vmapped.empty()) {
      return var(0.0);
    }

    // This code doesn't work with sliced arguments yet!
    const std::size_t num_vars_per_term = count_vars(vmapped[0]);
    const std::size_t num_vars_sliced_terms = num_terms * num_vars_per_term;
    const std::size_t num_vars_shared_terms = count_vars(args...);

    vari** varis = ChainableStack::instance_->memalloc_.alloc_array<vari*>(num_vars_sliced_terms + num_vars_shared_terms);
    double* partials = ChainableStack::instance_->memalloc_.alloc_array<double>(num_vars_sliced_terms + num_vars_shared_terms);

    save_varis(varis, vmapped);
    save_varis(varis + num_vars_sliced_terms, args...);

    /**
     * `local_stuff` is `tbb::enumerable_thread_specific<stuff>`
     *
     * On construction for each thread, build a ScopedChainableStack
     *  and make a copy of the args from the main stack
     *
     * Construction doesn't happen here, but inside `recursive_reducer`
     *  when the `.local` is accessed (I think that's when. At least
     *  it is not immediate).
     */
    local_stuff tlstuff([&]() {
      return stuff(new ScopedChainableStack(), args...);
    });

    recursive_reducer worker(tlstuff, std::forward<Vec>(vmapped), msgs);

    if (auto_partitioning) {
      tbb::parallel_reduce(
          tbb::blocked_range<std::size_t>(0, num_terms, grainsize), worker);
    } else {
      tbb::simple_partitioner partitioner;
      tbb::parallel_deterministic_reduce(
          tbb::blocked_range<std::size_t>(0, num_terms, grainsize), worker,
          partitioner);
    }

    /**
     * Count how many partial_sums there were. We will need to
     * save the vars for all these partial sums for the reverse
     * mode pass.
     */
    size_t M = 0;
    for (auto it = tlstuff.begin(); it != tlstuff.end(); ++it) {
      M += it->tails_.size();
    }

    ScopedChainableStack** stacks =
      ChainableStack::instance_->memalloc_.alloc_array<ScopedChainableStack*>(tlstuff.size());
    args_tuple_t** args_tuples =
      ChainableStack::instance_->memalloc_.alloc_array<args_tuple_t*>(tlstuff.size());
    var* tails = ChainableStack::instance_->memalloc_.alloc_array<var>(M);

    /**
     * Here we are copying all the information out of the
     * `stuff` data structures that we used in the forward pass
     * and are packing it up to be put in another data structure
     * for reverse pass.
     */
    size_t j = 0;
    size_t k = 0;
    for (auto it = tlstuff.begin(); it != tlstuff.end(); ++it) {
      stacks[j] = it->stack_; // per-thread ScopedChainableStack pointer
      args_tuples[j] = it->args_tuple_; // per-thread arg_tuple_ pointer
      for (size_t l = 0; l < it->tails_.size(); ++l) {
	tails[k] = it->tails_[l]; // per-partial-sum output adjoint
	k++;
      }
      j++;
    }

    /**
     * This is the vari attached to the `reduce_sum` operation that will live
     * on the main autodiff stack.
     *
     * It is responsible only for forwarding its adjoint along to another
     * data structure that actually does most of the work
     *
     * The reason this is broken into two pieces is because we need to have something
     * that has both a custom chain and a custom set_zero_all_adjoints function
     *
     * Because var must be constructed from a vari, we need a vari, and because
     * a vari cannot have a custom set_zero_all_adjoints, we need a second thing
     */
    struct reduce_sum_vari : public vari {
      reduce_sum_vari_impl* vb_;
      reduce_sum_vari(double val,
		      reduce_sum_vari_base* vb) : vari(val), vb_(vb) {}

      inline void chain() final {
	vb_->chain_(adj_);
      }
    };

    /**
     * This is the data structure that holds all the data to do the reverse
     * mode pass. It contains a lot of stuff that was used in the
     * `stuff` structure in the forward pass.
     *
     * The args_tuples_ argument needs a destructor called, which
     * is why this inherits from chainablle_alloc
     */
    struct reduce_sum_vari_base : public vari_base, public chainable_alloc {
      size_t N_;
      size_t M_;
      ScopedChainableStack** stacks_;
      args_tuple_t** args_tuples_;
      var* tails_;
      vari** varis_;
      size_t num_vars_shared_terms_;

      explicit reduce_sum_vari_base(size_t N, ScopedChainableStack** stacks, args_tuple_t** args_tuples, size_t M, var* tails, size_t num_vars_shared_terms, vari** varis)
	: N_(N), stacks_(stacks), args_tuples_(args_tuples),
	  M_(M), tails_(tails), num_vars_shared_terms_(num_vars_shared_terms),
	  varis_(varis) {
	ChainableStack::instance_->var_stack_.push_back(this);
      }

      ~reduce_sum_vari_base() {
	for(size_t i = 0; i < N_; ++i) {
	  delete stacks_[i];
	  delete args_tuples_[i];
	}
      }

      inline void chain() {}
      
      inline void chain_(double adj) {
	Eigen::VectorXd partials(num_vars_shared_terms_);
	partials.setZero();
	// Increment the tail of every term of the partial sum
	// with its adjoint
	for(size_t i = 0; i < M_; ++i) {
	  tails_[i].vi_->adj_ += adj;
	}
	for(size_t i = 0; i < N_; ++i) {
	  // Chain on every threaded stack -- doing this serially for laziness
	  stacks_[i]->execute([]() { grad(); });
	  // Pull the adjoints out and accumulate them in a temporary
	  apply([&](auto&&... args) {
	    accumulate_adjoints(partials.data(),
                                std::forward<decltype(args)>(args)...);
          }, *args_tuples_[i]);
	}
	// Accumulate the adjoints back on the main stack
	for(size_t j = 0; j < num_vars_shared_terms_; ++j) {
	  varis_[j]->adj_ += partials[j];
	}
      }

      // Setting all adjoints zero means on the threaded stacks too
      inline void set_zero_adjoint() final {
	for(size_t i = 0; i < N_; ++i) {
	  stacks_[i]->execute([]() { set_zero_all_adjoints(); });
	}
      }

      inline void init_dependent() {
      }
    };
    
    reduce_sum_vari_base* vb =
      new reduce_sum_vari_base(tlstuff.size(), stacks, args_tuples, M, tails, num_vars_shared_terms, varis);

    // This is the var that is the final result of our calculation!
    var res(new reduce_sum_vari(worker.sum_, vb));

    return res;
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
