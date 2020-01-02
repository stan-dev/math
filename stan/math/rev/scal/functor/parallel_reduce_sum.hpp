#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <iostream>
#include <iterator>
#include <vector>

namespace stan {
namespace math {

namespace internal {

template <class ReduceFunction, class M, class T, class Arg1, class Arg2>
struct parallel_reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2, var> {
  using vmapped_t = std::vector<M>;
  using arg1_t = std::vector<Arg1>;
  using arg2_t = std::vector<Arg2>;

  using vmapped_value_t = std::vector<typename child_type<M>::type>;
  using arg1_value_t = std::vector<typename child_type<Arg1>::type>;
  // using arg2_value_t = std::vector<typename child_type<Arg2>::type>;
  using arg2_value_t = std::vector<int>;

  using ops_partials_t
      = operands_and_partials<vmapped_t, std::vector<Arg1>, std::vector<Arg2>>;

  struct recursive_reducer {
    const vmapped_t& vmapped_;
    const arg1_t& arg1_;
    const arg1_value_t& arg1_value_;
    const arg2_t& arg2_;
    const arg2_value_t& arg2_value_;

    ops_partials_t& terms_partials_mapped_;
    ops_partials_t terms_partials_args_;
    double terms_sum_;

    recursive_reducer(const vmapped_t& vmapped, const T& init,
                      const arg1_t& arg1, const arg1_value_t& arg1_value,
                      const arg2_t& arg2, const arg2_value_t& arg2_value,
                      ops_partials_t& terms_partials_mapped)
        : vmapped_(vmapped),
          arg1_(arg1),
          arg1_value_(arg1_value),
          arg2_(arg2),
          arg2_value_(arg2_value),
          terms_partials_mapped_(terms_partials_mapped),
          terms_partials_args_(vmapped, arg1, arg2),
          terms_sum_(value_of(init)) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          arg1_(other.arg1_),
          arg1_value_(other.arg1_value_),
          arg2_(other.arg2_),
          arg2_value_(other.arg2_value_),
          terms_partials_mapped_(other.terms_partials_mapped_),
          terms_partials_args_(vmapped_, arg1_, arg2_),
          terms_sum_(0.0) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      try {
        start_nested();

        // elem_t should better be a non-var!!! Otherwise we interfere
        // with the outer AD tree as we mess with the adjoints (and
        // right now we discard that it is a var if it would be)
        vmapped_t local_sub_slice;

        local_sub_slice.reserve(r.size());
        for (std::size_t i = r.begin(); i != r.end(); ++i)
          local_sub_slice.emplace_back(value_of(vmapped_[i]));

        // create a deep copy of arg1 which is not tied to the outer
        // AD tree
        // todo: these copies can be put onto a thread-local storage
        // object such that we only create these once per thread which
        // should significantly boost performance as many copies are
        // avoided... this comes at the cost that the adjoints have to
        // be zeroed "manually"
        /* using the var_nochain_stack_ tape avoids up-propagation of chain
        std::vector<Arg1> varg1;
        varg1.reserve(varg1_.size());

        for (Arg1& elem : varg1_)
          varg1.emplace_back(var(new vari(elem.val(), false)));
        */

        // but this should actually do the same if we simply
        // instantiate the var from a double representation
        // const std::vector<double> varg1_d = value_of(varg1_);
        // std::vector<Arg1> varg1(varg1_d.begin(), varg1_d.end());

        arg1_t local_arg1(arg1_value_.begin(), arg1_value_.end());
        arg2_t local_arg2(arg2_value_.begin(), arg2_value_.end());

        T sub_sum_v = ReduceFunction()(r.begin(), r.end() - 1, local_sub_slice,
                                       local_arg1, local_arg2);

        sub_sum_v.grad();

        terms_sum_ += sub_sum_v.val();

        /*
        if (!is_constant_all<M>::value) {
          for (std::size_t i = 0; i != r.size(); ++i)
            terms_partials_mapped_.edge1_.partials_[r.begin() + i] +=
        local_sub_slice[i].adj();
        }
        */

        if (!is_constant_all<Arg1>::value) {
          for (std::size_t i = 0; i != local_arg1.size(); ++i)
            terms_partials_args_.edge2_.partials_[i] += local_arg1[i].adj();
        }

        /*
        if (!is_constant_all<Arg2>::value) {
          for (std::size_t i = 0; i != local_arg2.size(); ++i)
            terms_partials_args_.edge3_.partials_[i] += local_arg2[i].adj();
        }
        */

      } catch (const std::exception& e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
    }

    void join(const recursive_reducer& child) {
      terms_sum_ += child.terms_sum_;

      if (!is_constant_all<Arg1>::value) {
        for (std::size_t i = 0; i != arg1_.size(); ++i)
          terms_partials_args_.edge2_.partials_[i]
              += child.terms_partials_args_.edge2_.partials_[i];
      }

      /*
      if (!is_constant_all<Arg2>::value) {
        for (std::size_t i = 0; i != arg2_.size(); ++i)
          terms_partials_args_.edge3_.partials_[i]
              += child.terms_partials_args_.edge3_.partials_[i];
      }
      */
    }
  };

  T operator()(const vmapped_t& vmapped, T init, std::size_t grainsize,
               const arg1_t& arg1, const arg2_t& arg2) const {
    const std::size_t num_jobs = vmapped.size();
    arg1_value_t arg1_value = value_of(arg1);
    arg2_value_t arg2_value = value_of(arg2);

    ops_partials_t ops(vmapped, arg1, arg2);

    recursive_reducer worker(vmapped, init, arg1, arg1_value, arg2, arg2_value,
                             ops);

    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);

    if (!is_constant_all<Arg1>::value) {
      ops.edge2_.partials_ = worker.terms_partials_args_.edge2_.partials_;
    }

    if (!is_constant_all<Arg2>::value) {
      ops.edge3_.partials_ = worker.terms_partials_args_.edge3_.partials_;
    }

    return ops.build(worker.terms_sum_);
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
