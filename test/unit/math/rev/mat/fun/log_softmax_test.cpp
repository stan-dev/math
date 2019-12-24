#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

TEST(AgradRevMatrix, log_softmax_nested) {
  using stan::math::log_softmax;
  using stan::math::row_vector_v;
  using stan::math::vector_v;
  using stan::math::vector_d;
  using stan::math::var;

  std::vector<vector_v> v(5);
  std::vector<row_vector_v> rv(5);
  std::vector<std::vector<var>> stv(5);

  std::vector<vector_v> v_result_iter(5);
  std::vector<row_vector_v> rv_result_iter(5);
  std::vector<std::vector<var>> stv_result_iter(5);

  for(int i = 0; i < 5; ++i){
    v[i] = vector_d::Random(10);
    v[i].adj() = vector_d::Random(10);
    rv[i] = v[i].transpose();
    stv[i] = std::vector<var>(v[i].data(), v[i].data() + v[i].size());

    v_result_iter[i] = log_softmax(v[i]);
    rv_result_iter[i] = log_softmax(rv[i]);
    stv_result_iter[i] = log_softmax(stv[i]);
  }

  std::vector<vector_v> v_result = log_softmax(v);
  std::vector<row_vector_v> rv_result = log_softmax(rv);
  std::vector<std::vector<var>> stv_result = log_softmax(stv);

  for(int i = 0; i < 5; ++i){
    Eigen::Map<vector_v> stv_result_map(stv_result[i].data(),
                                        stv_result[i].size());
    Eigen::Map<vector_v> stv_result_iter_map(stv_result_iter[i].data(),
                                             stv_result_iter[i].size());
    expect_matrix_eq(v_result[i].val().array() * 3,
                     v_result_iter[i].val() +
                        rv_result_iter[i].val().transpose() + 
                        stv_result_iter_map.val());
    expect_matrix_eq(v_result[i].adj().array() * 3,
                     v_result_iter[i].adj() +
                        rv_result_iter[i].adj().transpose() + 
                        stv_result_iter_map.adj());
    expect_matrix_eq(v_result[i].val().array() * 2,
                     rv_result[i].val().transpose() + 
                        stv_result_map.val());
    expect_matrix_eq(v_result[i].adj().array() * 2,
                     rv_result[i].adj().transpose() + 
                        stv_result_map.adj());
  }
}
