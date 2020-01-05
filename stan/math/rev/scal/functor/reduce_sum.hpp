#ifndef STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_REDUCE_SUM_HPP

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

template <>
struct child_type<int> {
  using type = int;
};

template <typename T>
typename child_type<T>::type adjoint_of(const T x) {
  return typename child_type<T>::type(0);
}

template <typename T>
std::vector<typename child_type<T>::type> adjoint_of(const std::vector<T>& x) {
  const size_t size = x.size();
  std::vector<typename child_type<T>::type> result(size);
  for (size_t i = 0; i != size; ++i) {
    result[i] = adjoint_of(x[i]);
  }
  return result;
}

double adjoint_of(const var& x) { return x.vi_->adj_; }

namespace internal {

template <class ReduceFunction, class M, class T, class Arg1, class Arg2,
          class Arg3, class Arg4>
struct reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2, Arg3, Arg4, var> {
  using vmapped_t = std::vector<M>;
  using arg1_t = std::vector<Arg1>;
  using arg2_t = std::vector<Arg2>;
  using arg3_t = std::vector<Arg3>;
  using arg4_t = std::vector<Arg4>;

  using vmapped_value_t = std::vector<typename child_type<M>::type>;
  using arg1_value_t = std::vector<typename child_type<Arg1>::type>;
  using arg2_value_t = std::vector<typename child_type<Arg2>::type>;
  using arg3_value_t = std::vector<typename child_type<Arg3>::type>;
  using arg4_value_t = std::vector<typename child_type<Arg4>::type>;

  using ops_partials_t
      = operands_and_partials<vmapped_t, std::vector<Arg1>, std::vector<Arg2>,
                              std::vector<Arg3>, std::vector<Arg4>>;

  struct recursive_reducer {
    const vmapped_t& vmapped_;
    const arg1_t& arg1_;
    const arg1_value_t& arg1_value_;
    const arg2_t& arg2_;
    const arg2_value_t& arg2_value_;
    const arg3_t& arg3_;
    const arg3_value_t& arg3_value_;
    const arg4_t& arg4_;
    const arg4_value_t& arg4_value_;

    ops_partials_t& terms_partials_mapped_;
    ops_partials_t terms_partials_args_;
    double terms_sum_;

    recursive_reducer(const vmapped_t& vmapped, const T& init,
                      const arg1_t& arg1, const arg1_value_t& arg1_value,
                      const arg2_t& arg2, const arg2_value_t& arg2_value,
                      const arg3_t& arg3, const arg3_value_t& arg3_value,
                      const arg4_t& arg4, const arg4_value_t& arg4_value,
                      ops_partials_t& terms_partials_mapped)
        : vmapped_(vmapped),
          arg1_(arg1),
          arg1_value_(arg1_value),
          arg2_(arg2),
          arg2_value_(arg2_value),
          arg3_(arg3),
          arg3_value_(arg3_value),
          arg4_(arg4),
          arg4_value_(arg4_value),
          terms_partials_mapped_(terms_partials_mapped),
          terms_partials_args_(vmapped, arg1, arg2, arg3, arg4),
          terms_sum_(value_of(init)) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          arg1_(other.arg1_),
          arg1_value_(other.arg1_value_),
          arg2_(other.arg2_),
          arg2_value_(other.arg2_value_),
          arg3_(other.arg3_),
          arg3_value_(other.arg3_value_),
          arg4_(other.arg4_),
          arg4_value_(other.arg4_value_),
          terms_partials_mapped_(other.terms_partials_mapped_),
          terms_partials_args_(vmapped_, arg1_, arg2_, arg3_, arg4_),
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

        const arg1_t local_arg1(arg1_value_.begin(), arg1_value_.end());
        const arg2_t local_arg2(arg2_value_.begin(), arg2_value_.end());
        const arg3_t local_arg3(arg3_value_.begin(), arg3_value_.end());
        const arg4_t local_arg4(arg4_value_.begin(), arg4_value_.end());

        T sub_sum_v
            = ReduceFunction()(r.begin(), r.end() - 1, local_sub_slice,
                               local_arg1, local_arg2, local_arg3, local_arg4);

        sub_sum_v.grad();

        terms_sum_ += sub_sum_v.val();

        if (!is_constant_all<M>::value) {
          vmapped_value_t sub_slice_adjoint = adjoint_of(local_sub_slice);
          for (std::size_t i = 0; i != r.size(); ++i)
            terms_partials_mapped_.edge1_.partials_[r.begin() + i]
                += sub_slice_adjoint[i];
        }
        if (!is_constant_all<Arg1>::value) {
          arg1_value_t arg1_adjoint = adjoint_of(local_arg1);
          for (std::size_t i = 0; i != local_arg1.size(); ++i)
            terms_partials_args_.edge2_.partials_[i] += arg1_adjoint[i];
        }
        if (!is_constant_all<Arg2>::value) {
          arg2_value_t arg2_adjoint = adjoint_of(local_arg2);
          for (std::size_t i = 0; i != local_arg2.size(); ++i)
            terms_partials_args_.edge3_.partials_[i] += arg2_adjoint[i];
        }
        if (!is_constant_all<Arg3>::value) {
          arg3_value_t arg3_adjoint = adjoint_of(local_arg3);
          for (std::size_t i = 0; i != local_arg3.size(); ++i)
            terms_partials_args_.edge4_.partials_[i] += arg3_adjoint[i];
        }
        if (!is_constant_all<Arg4>::value) {
          arg4_value_t arg4_adjoint = adjoint_of(local_arg4);
          for (std::size_t i = 0; i != local_arg2.size(); ++i)
            terms_partials_args_.edge5_.partials_[i] += arg4_adjoint[i];
        }

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
      if (!is_constant_all<Arg2>::value) {
        for (std::size_t i = 0; i != arg2_.size(); ++i)
          terms_partials_args_.edge3_.partials_[i]
              += child.terms_partials_args_.edge3_.partials_[i];
      }
      if (!is_constant_all<Arg3>::value) {
        for (std::size_t i = 0; i != arg3_.size(); ++i)
          terms_partials_args_.edge4_.partials_[i]
              += child.terms_partials_args_.edge4_.partials_[i];
      }
      if (!is_constant_all<Arg4>::value) {
        for (std::size_t i = 0; i != arg4_.size(); ++i)
          terms_partials_args_.edge5_.partials_[i]
              += child.terms_partials_args_.edge5_.partials_[i];
      }
    }
  };

  T operator()(const vmapped_t& vmapped, T init, std::size_t grainsize,
               const arg1_t& arg1, const arg2_t& arg2, const arg3_t& arg3,
               const arg4_t& arg4) const {
    const std::size_t num_jobs = vmapped.size();
    arg1_value_t arg1_value = value_of(arg1);
    arg2_value_t arg2_value = value_of(arg2);
    arg3_value_t arg3_value = value_of(arg3);
    arg4_value_t arg4_value = value_of(arg4);

    ops_partials_t ops(vmapped, arg1, arg2, arg3, arg4);

    recursive_reducer worker(vmapped, init, arg1, arg1_value, arg2, arg2_value,
                             arg3, arg3_value, arg4, arg4_value, ops);

    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);

    if (!is_constant_all<Arg1>::value) {
      ops.edge2_.partials_ = worker.terms_partials_args_.edge2_.partials_;
    }
    if (!is_constant_all<Arg2>::value) {
      ops.edge3_.partials_ = worker.terms_partials_args_.edge3_.partials_;
    }
    if (!is_constant_all<Arg3>::value) {
      ops.edge4_.partials_ = worker.terms_partials_args_.edge4_.partials_;
    }
    if (!is_constant_all<Arg4>::value) {
      ops.edge5_.partials_ = worker.terms_partials_args_.edge5_.partials_;
    }

    return ops.build(worker.terms_sum_);
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
