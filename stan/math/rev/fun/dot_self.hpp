#ifndef STAN_MATH_REV_FUN_DOT_SELF_HPP
#define STAN_MATH_REV_FUN_DOT_SELF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <vector>

namespace stan {
namespace math {

namespace internal {	
class dot_self_vari : public vari {	
protected:	
  vari** v_;	
  size_t size_;	

public:	
  dot_self_vari(vari** v, size_t size)	
      : vari(Eigen::Map<vector_vi>(v, size).val().squaredNorm()),	
        v_(v),	
        size_(size) {}	
  template <typename T, require_eigen_t<T>* = nullptr>	
  explicit dot_self_vari(const T& v)	
      : vari(v.val().squaredNorm()), size_(v.size()) {	
    v_ = reinterpret_cast<vari**>(	
        ChainableStack::instance_->memalloc_.alloc(size_ * sizeof(vari*)));	
    Eigen::Map<matrix_vi>(v_, v.rows(), v.cols()) = v.vi();	
  }	
  virtual void chain() {	
    Eigen::Map<vector_vi> v_map(v_, size_);	
    v_map.adj() += adj_ * 2.0 * v_map.val();	
  }	
};	
}
  
/**
 * Returns the dot product of a vector of var with itself.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase and
 * have one compile time dimension equal to 1)
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename T, require_eigen_vector_vt<is_var, T>* = nullptr>
inline var dot_self(const T& v) {
  /*
    --------------------------------------------------------------------
    Benchmark                          Time             CPU   Iterations
    --------------------------------------------------------------------
    toss_me                     10126264 ns     10125955 ns           70
    dot_self/2/manual_time          65.2 ns          161 ns     10681302
    dot_self/4/manual_time          66.7 ns          190 ns     10469218
    dot_self/8/manual_time          71.7 ns          260 ns      9693894
    dot_self/16/manual_time         84.3 ns          378 ns      8305816
    dot_self/32/manual_time          121 ns          635 ns      5814113
    dot_self/64/manual_time          186 ns         1119 ns      3757719
    dot_self/128/manual_time         326 ns         2105 ns      2150730
    dot_self/256/manual_time         611 ns         4212 ns      1142432
    dot_self/512/manual_time        1217 ns         8430 ns       574490
    dot_self/1024/manual_time       2548 ns        19862 ns       275469
  */

  //const Eigen::Ref<const plain_type_t<T>>& v_ref = v;
  //return {new internal::dot_self_vari(v_ref)};
  
  /*
    --------------------------------------------------------------------
    Benchmark                          Time             CPU   Iterations
    --------------------------------------------------------------------
    toss_me                     10077101 ns     10076684 ns           69
    dot_self/2/manual_time          68.7 ns          164 ns     10057031
    dot_self/4/manual_time          69.2 ns          192 ns     10049326
    dot_self/8/manual_time          74.2 ns          263 ns      9370878
    dot_self/16/manual_time         89.4 ns          385 ns      7810417
    dot_self/32/manual_time          125 ns          671 ns      5615422
    dot_self/64/manual_time          195 ns         1149 ns      3584751
    dot_self/128/manual_time         348 ns         2166 ns      2008505
    dot_self/256/manual_time         651 ns         4508 ns      1069752
    dot_self/512/manual_time        1277 ns         8948 ns       544405
    dot_self/1024/manual_time       2648 ns        21018 ns       267065
  */
  
  arena_t<T> arena_v = v;

  var res = arena_v.val().squaredNorm();

  reverse_pass_callback([arena_v, res]() {
    for(size_t i = 0; i < arena_v.size(); ++i) {
      arena_v.coeffRef(i).adj() += 2.0 * res.adj() * arena_v.coeff(i).val();
    }
  });

  return res;
  

  /*
    --------------------------------------------------------------------
    Benchmark                          Time             CPU   Iterations
    --------------------------------------------------------------------
    toss_me                     10161766 ns     10160959 ns           70
    dot_self/2/manual_time          84.7 ns          180 ns      8201199
    dot_self/4/manual_time          85.6 ns          210 ns      8076649
    dot_self/8/manual_time          92.1 ns          286 ns      7451948
    dot_self/16/manual_time          103 ns          422 ns      6857630
    dot_self/32/manual_time          123 ns          675 ns      5680814
    dot_self/64/manual_time          172 ns         1210 ns      4071015
    dot_self/128/manual_time         276 ns         2272 ns      2544438
    dot_self/256/manual_time         505 ns         4629 ns      1384535
    dot_self/512/manual_time         992 ns         9303 ns       703375
    dot_self/1024/manual_time       2270 ns        20923 ns       308428
  */
  
  /*
  arena_t<T> arena_v = v;

  auto arena_v_val = to_arena(value_of(arena_v));
  var res = arena_v_val.squaredNorm();
  
  reverse_pass_callback([arena_v, res, arena_v_val]() mutable {
  arena_v.adj() += (2.0 * res.adj()) * arena_v_val;
  });
  
  return res;
  */
}

}  // namespace math
}  // namespace stan
#endif
