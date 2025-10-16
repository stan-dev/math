#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>
#include <limits>
#include <type_traits>
#include <vector>

template <typename T>
auto one_arg_bad_vals(const T& x) {
  return x;
}

template <typename T>
auto one_arg_bad_vals(const stan::math::var_value<T>& x) {
  return stan::math::make_callback_var(
      x.val() * 0, [x](const auto& vi) mutable { x.adj() += vi.adj(); });
}

TEST(test_unit_math_test_ad_matvar, one_arg_bad_vals) {
  auto f = [](const auto& u) { return one_arg_bad_vals(u); };

  Eigen::VectorXd x = Eigen::VectorXd::Ones(2);

  // There are two failures here and Google test requires
  //  that we expect them both:
  //  https://stackoverflow.com/questions/24037390/how-can-i-expect-multiple-failures-in-google-test
  //  https://github.com/google/googletest/blob/35fb11efbe1a2761ce923f49a9df1a430e5d16be/googletest/docs/AdvancedGuide.md#catching-failures
  // They do not work properly unless you also include a string to match (the
  // outer expect will ignore the inner expect's failure to fail)
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x), "values"),
      "values");
}

template <typename T>
auto one_arg_bad_grads(const T& x) {
  return x;
}

template <typename T>
auto one_arg_bad_grads(const stan::math::var_value<T>& x) {
  return stan::math::make_callback_var(
      x.val(), [x](const auto& vi) mutable { x.adj() -= vi.adj(); });
}

TEST(test_unit_math_test_ad_matvar, one_arg_bad_grads) {
  auto f = [](const auto& u) { return one_arg_bad_grads(u); };

  Eigen::VectorXd x = Eigen::VectorXd::Ones(2);

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x), "adjoints"),
      "adjoints");
}

template <typename T1, typename T2>
auto two_arg_bad_vals(const T1& x1, const T2& x2) {
  return stan::math::add(x1, x2);
}

template <typename T1, typename T2, stan::require_st_arithmetic<T1>* = nullptr>
auto two_arg_bad_vals(const T1& x1, const stan::math::var_value<T2>& x2) {
  Eigen::VectorXd ret_val = 0.0 * stan::math::add(x1, x2.val());
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1,
                                             stan::math::var_value<T2>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x2, ret]() mutable {
    if (stan::is_stan_scalar<T2>::value) {
      stan::math::forward_as<stan::math::var>(x2).adj()
          += stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x2.adj()) += ret.adj();
    }
  });

  return ret_type(ret);
}

template <typename T1, typename T2, stan::require_st_arithmetic<T1>* = nullptr>
auto two_arg_bad_vals(const stan::math::var_value<T2>& x2, const T1& x1) {
  Eigen::VectorXd ret_val = 0.0 * stan::math::add(x1, x2.val());
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1,
                                             stan::math::var_value<T2>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x2, ret]() mutable {
    if (stan::is_stan_scalar<T2>::value) {
      stan::math::forward_as<stan::math::var>(x2).adj()
          += stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x2.adj()) += ret.adj();
    }
  });

  return ret_type(ret);
}

TEST(test_unit_math_test_ad_matvar, two_arg_bad_vals) {
  auto f
      = [](const auto& x1, const auto& x2) { return two_arg_bad_vals(x1, x2); };

  Eigen::VectorXd x = Eigen::VectorXd::Ones(2);
  double y = 1.0;

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, x), "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, y), "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, x), "values"),
      "values");
}

template <typename T1, typename T2>
auto two_arg_bad_grads(const T1& x1, const T2& x2) {
  return stan::math::add(x1, x2);
}

template <typename T1, typename T2, stan::require_st_arithmetic<T1>* = nullptr>
auto two_arg_bad_grads(const T1& x1, const stan::math::var_value<T2>& x2) {
  auto ret_val = stan::math::add(x1, x2.val());
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1,
                                             stan::math::var_value<T2>>;
  stan::arena_t<ret_type> ret = stan::math::add(x1, x2.val());

  stan::math::reverse_pass_callback([x2, ret]() mutable {
    if (stan::is_stan_scalar<T2>::value) {
      stan::math::forward_as<stan::math::var>(x2).adj()
          -= stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x2.adj()) -= ret.adj();
    }
  });

  return ret_type(ret);
}

template <typename T1, typename T2, stan::require_st_arithmetic<T1>* = nullptr>
auto two_arg_bad_grads(const stan::math::var_value<T2>& x2, const T1& x1) {
  auto ret_val = stan::math::add(x1, x2.val());
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1,
                                             stan::math::var_value<T2>>;
  stan::arena_t<ret_type> ret = stan::math::add(x1, x2.val());

  stan::math::reverse_pass_callback([x2, ret]() mutable {
    if (stan::is_stan_scalar<T2>::value) {
      stan::math::forward_as<stan::math::var>(x2).adj()
          -= stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x2.adj()) -= ret.adj();
    }
  });

  return ret_type(ret);
}

TEST(test_unit_math_test_ad_matvar, two_arg_bad_grads) {
  auto f = [](const auto& x1, const auto& x2) {
    return two_arg_bad_grads(x1, x2);
  };

  Eigen::VectorXd x = Eigen::VectorXd::Ones(2);
  double y = 1.0;

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, x),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, y),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, x),
                              "adjoints"),
      "adjoints");
}

template <typename T1, typename T2, typename T3>
auto three_arg_bad_vals(const T1& x1, const T2& x2, const T3& x3) {
  return stan::math::add(x1, stan::math::add(x2, x3));
}

template <typename T1, typename T2, typename T3,
          stan::require_all_st_arithmetic<T1, T2>* = nullptr>
auto three_arg_bad_vals(const T1& x1, const T2& x2,
                        const stan::math::var_value<T3>& x3) {
  auto ret_val = stan::math::add(x1, stan::math::add(x2, 0.0 * x3.val()));
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1, T2,
                                             stan::math::var_value<T3>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x3, ret]() mutable {
    if (stan::is_stan_scalar<T3>::value) {
      stan::math::forward_as<stan::math::var>(x3).adj()
          += stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x3.adj()) += ret.adj();
    }
  });

  return ret_type(ret);
}

template <typename T1, typename T2, typename T3,
          stan::require_all_st_arithmetic<T1, T2>* = nullptr>
auto three_arg_bad_vals(const T1& x1, const stan::math::var_value<T3>& x3,
                        const T2& x2) {
  auto ret_val = stan::math::add(x1, stan::math::add(x2, 0.0 * x3.val()));
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1, T2,
                                             stan::math::var_value<T3>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x3, ret]() mutable {
    if (stan::is_stan_scalar<T3>::value) {
      stan::math::forward_as<stan::math::var>(x3).adj()
          += stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x3.adj()) += ret.adj();
    }
  });

  return ret_type(ret);
}

template <typename T1, typename T2, typename T3,
          stan::require_all_st_arithmetic<T1, T2>* = nullptr>
auto three_arg_bad_vals(const stan::math::var_value<T3>& x3, const T1& x1,
                        const T2& x2) {
  auto ret_val = stan::math::add(x1, stan::math::add(x2, 0.0 * x3.val()));
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1, T2,
                                             stan::math::var_value<T3>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x3, ret]() mutable {
    if (stan::is_stan_scalar<T3>::value) {
      stan::math::forward_as<stan::math::var>(x3).adj()
          += stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x3.adj()) += ret.adj();
    }
  });

  return ret_type(ret);
}

TEST(test_unit_math_test_ad_matvar, three_arg_bad_vals) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return three_arg_bad_vals(x1, x2, x3);
  };

  Eigen::VectorXd x = Eigen::VectorXd::Ones(2);
  double y = 1.0;

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, y, x),
                              "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, x, y),
                              "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, x, x),
                              "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, y, y),
                              "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, y, x),
                              "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, x, y),
                              "values"),
      "values");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, x, x),
                              "values"),
      "values");
}

template <typename T1, typename T2, typename T3>
auto three_arg_bad_grads(const T1& x1, const T2& x2, const T3& x3) {
  return stan::math::add(x1, stan::math::add(x2, x3));
}

template <typename T1, typename T2, typename T3,
          stan::require_all_st_arithmetic<T1, T2>* = nullptr>
auto three_arg_bad_grads(const T1& x1, const T2& x2,
                         const stan::math::var_value<T3>& x3) {
  auto ret_val = stan::math::add(x1, stan::math::add(x2, x3.val()));
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1, T2,
                                             stan::math::var_value<T3>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x3, ret]() mutable {
    if (stan::is_stan_scalar<T3>::value) {
      stan::math::forward_as<stan::math::var>(x3).adj()
          -= stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x3.adj()) -= ret.adj();
    }
  });

  return ret_type(ret);
}

template <typename T1, typename T2, typename T3,
          stan::require_all_st_arithmetic<T1, T2>* = nullptr>
auto three_arg_bad_grads(const T1& x1, const stan::math::var_value<T3>& x3,
                         const T2& x2) {
  auto ret_val = stan::math::add(x1, stan::math::add(x2, x3.val()));
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1, T2,
                                             stan::math::var_value<T3>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x3, ret]() mutable {
    if (stan::is_stan_scalar<T3>::value) {
      stan::math::forward_as<stan::math::var>(x3).adj()
          -= stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x3.adj()) -= ret.adj();
    }
  });

  return ret_type(ret);
}

template <typename T1, typename T2, typename T3,
          stan::require_all_st_arithmetic<T1, T2>* = nullptr>
auto three_arg_bad_grads(const stan::math::var_value<T3>& x3, const T1& x1,
                         const T2& x2) {
  auto ret_val = stan::math::add(x1, stan::math::add(x2, x3.val()));
  using ret_type = stan::return_var_matrix_t<decltype(ret_val), T1, T2,
                                             stan::math::var_value<T3>>;
  stan::arena_t<ret_type> ret = ret_val;

  stan::math::reverse_pass_callback([x3, ret]() mutable {
    if (stan::is_stan_scalar<T3>::value) {
      stan::math::forward_as<stan::math::var>(x3).adj()
          -= stan::math::sum(ret.adj());
    } else {
      stan::math::forward_as<Eigen::VectorXd>(x3.adj()) -= ret.adj();
    }
  });

  return ret_type(ret);
}

TEST(test_unit_math_test_ad_matvar, three_arg_bad_grads) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return three_arg_bad_grads(x1, x2, x3);
  };

  Eigen::VectorXd x = Eigen::VectorXd::Ones(2);
  double y = 1.0;

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, y, x),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, x, y),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, y, x, x),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, y, y),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, y, x),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, x, y),
                              "adjoints"),
      "adjoints");
  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x, x, x),
                              "adjoints"),
      "adjoints");
}

template <typename T>
auto one_arg_bad_vals_std_vector(const std::vector<T>& x) {
  return x;
}

template <typename T>
auto one_arg_bad_vals_std_vector(
    const std::vector<stan::math::var_value<T>>& x) {
  using ret_type = std::vector<stan::math::var_value<T>>;
  stan::arena_t<ret_type> arena_x(x.size());

  stan::arena_t<ret_type> out(x.size());
  ret_type out_heap(out.size());
  for (size_t i = 0; i < x.size(); ++i) {
    arena_x[i] = x[i];
    out[i] = x[i].val() * 0.0;
    out_heap[i] = out[i];
  }

  stan::math::reverse_pass_callback([arena_x, out]() mutable {
    for (size_t i = 0; i < arena_x.size(); ++i) {
      arena_x[i].adj() += out[i].adj();
    }
  });

  return out_heap;
}

TEST(test_unit_math_test_ad_matvar, one_arg_bad_vals_std_vector) {
  auto f = [](const auto& u) { return one_arg_bad_vals_std_vector(u); };

  std::vector<Eigen::VectorXd> x = {Eigen::VectorXd::Ones(2)};

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x), "values"),
      "values");
}

template <typename T>
auto one_arg_bad_grads_std_vector(const std::vector<T>& x) {
  return x;
}

template <typename T>
auto one_arg_bad_grads_std_vector(
    const std::vector<stan::math::var_value<T>>& x) {
  using ret_type = std::vector<stan::math::var_value<T>>;
  stan::arena_t<ret_type> arena_x(x.size());

  stan::arena_t<ret_type> out(x.size());
  ret_type out_heap(out.size());
  for (size_t i = 0; i < x.size(); ++i) {
    arena_x[i] = x[i];
    out[i] = x[i].val();
    out_heap[i] = out[i];
  }

  stan::math::reverse_pass_callback([arena_x, out]() mutable {
    for (size_t i = 0; i < arena_x.size(); ++i) {
      arena_x[i].adj() -= out[i].adj();
    }
  });

  return out_heap;
}

TEST(test_unit_math_test_ad_matvar, one_arg_bad_grads_std_vector) {
  auto f = [](const auto& u) { return one_arg_bad_grads_std_vector(u); };

  std::vector<Eigen::VectorXd> x = {Eigen::VectorXd::Ones(2)};

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x), "adjoints"),
      "adjoints");
}

template <typename T1, typename T2>
auto two_args_bad_vals_std_vector(const std::vector<T1>& x1,
                                  const std::vector<T2>& x2) {
  std::vector<stan::plain_type_t<decltype(
      stan::math::add(std::declval<T1>(), std::declval<T2>()))>>
      array_sum;
  for (size_t i = 0; i < x1.size(); ++i) {
    array_sum.push_back(stan::math::add(x1[i], x2[i]));
  }

  return array_sum;
}

template <typename T1, typename T2>
auto two_args_bad_vals_std_vector(
    const std::vector<T1>& x1,
    const std::vector<stan::math::var_value<T2>>& x2) {
  using ret_type = std::vector<stan::math::var_value<T2>>;
  stan::arena_t<std::vector<T1>> arena_x1(x1.size());
  stan::arena_t<ret_type> arena_x2(x2.size());

  stan::arena_t<ret_type> out(x1.size());
  ret_type out_heap(out.size());
  for (size_t i = 0; i < x1.size(); ++i) {
    arena_x1[i] = x1[i];
    arena_x2[i] = x2[i];
    out[i] = x2[i].val() * 0.0;
    out_heap[i] = out[i];
  }

  stan::math::reverse_pass_callback([arena_x1, arena_x2, out]() mutable {
    for (size_t i = 0; i < arena_x1.size(); ++i) {
      if (!stan::is_constant<T1>::value) {
        stan::math::forward_as<
            stan::return_var_matrix_t<T2, T1, stan::math::var>>(arena_x1[i])
            .adj()
            += out[i].adj();
      }
      arena_x2[i].adj() += out[i].adj();
    }
  });

  return out_heap;
}

TEST(test_unit_math_test_ad_matvar, two_args_bad_vals_std_vector) {
  auto f = [](const auto& x1, const auto& x2) {
    return two_args_bad_vals_std_vector(x1, x2);
  };

  std::vector<Eigen::VectorXd> x1 = {Eigen::VectorXd::Ones(2)};
  std::vector<Eigen::VectorXd> x2 = {Eigen::VectorXd::Zero(2)};

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x1, x2);
                              , "values"),
      "values");
}

template <typename T1, typename T2>
auto two_args_bad_grads_std_vector(const std::vector<T1>& x1,
                                   const std::vector<T2>& x2) {
  std::vector<stan::plain_type_t<decltype(
      stan::math::add(std::declval<T1>(), std::declval<T2>()))>>
      array_sum;
  for (size_t i = 0; i < x1.size(); ++i) {
    array_sum.push_back(stan::math::add(x1[i], x2[i]));
  }

  return array_sum;
}

template <typename T1, typename T2>
auto two_args_bad_grads_std_vector(
    const std::vector<T1>& x1,
    const std::vector<stan::math::var_value<T2>>& x2) {
  using ret_type = std::vector<stan::math::var_value<T2>>;
  stan::arena_t<std::vector<T1>> arena_x1(x1.size());
  stan::arena_t<ret_type> arena_x2(x2.size());

  stan::arena_t<ret_type> out(x1.size());
  ret_type out_heap(out.size());
  for (size_t i = 0; i < x1.size(); ++i) {
    arena_x1[i] = x1[i];
    arena_x2[i] = x2[i];
    out[i] = stan::math::add(stan::math::value_of(x1[i]), x2[i].val());
    out_heap[i] = out[i];
  }

  stan::math::reverse_pass_callback([arena_x1, arena_x2, out]() mutable {
    for (size_t i = 0; i < arena_x1.size(); ++i) {
      if (!stan::is_constant<T1>::value) {
        stan::math::forward_as<
            stan::return_var_matrix_t<T2, T1, stan::math::var>>(arena_x1[i])
            .adj()
            += out[i].adj();
      }
      arena_x2[i].adj() += out[i].adj() * 0;
    }
  });

  return out_heap;
}

TEST(test_unit_math_test_ad_matvar, two_args_bad_grads_std_vector) {
  auto f = [](const auto& x1, const auto& x2) {
    return two_args_bad_grads_std_vector(x1, x2);
  };

  std::vector<Eigen::VectorXd> x1 = {Eigen::VectorXd::Ones(2)};
  std::vector<Eigen::VectorXd> x2 = {Eigen::VectorXd::Zero(2)};

  EXPECT_NONFATAL_FAILURE(
      EXPECT_NONFATAL_FAILURE(stan::test::expect_ad_matvar(f, x1, x2);
                              , "adjoints"),
      "adjoints");
}

template <typename T>
auto bad_return_type(const T& x) {
  return x;
}

template <typename T>
auto bad_return_type(const stan::math::var_value<T>& x) {
  using ret_type = stan::math::promote_scalar_t<stan::math::var, T>;
  stan::arena_t<ret_type> out = x.val();

  stan::math::reverse_pass_callback(
      [x, out]() mutable { x.adj() += out.adj(); });

  return ret_type(out);
}

TEST(test_unit_math_test_ad_matvar, bad_return_type) {
  EXPECT_FATAL_FAILURE(
      {
        auto f = [](const auto& x) { return bad_return_type(x); };

        Eigen::VectorXd x = Eigen::VectorXd::Ones(2);

        stan::test::expect_ad_matvar(f, x);
      },
      "");
}

template <typename T>
auto bad_return_type_std_vector(const std::vector<T>& x) {
  return x;
}

template <typename T>
auto bad_return_type_std_vector(
    const std::vector<stan::math::var_value<T>>& x) {
  std::vector<stan::math::promote_scalar_t<stan::math::var, T>> out(x.size());

  for (size_t i = 0; i < out.size(); ++i) {
    out[i] = x[i].val();
  }

  return out;
}

TEST(test_unit_math_test_ad_matvar, bad_return_type_std_vector) {
  EXPECT_FATAL_FAILURE(
      {
        auto f = [](const auto& x) { return bad_return_type_std_vector(x); };

        std::vector<Eigen::VectorXd> x = {Eigen::VectorXd::Ones(2)};

        stan::test::expect_ad_matvar(f, x);
      },
      "");
}

template <typename T>
auto bad_throws1(const T& x) {
  return x;
}

template <typename T>
auto bad_throws1(const stan::math::var_value<T>& x) {
  throw std::exception();
  return x;
}

TEST(test_unit_math_test_ad_matvar, bad_throws1) {
  EXPECT_FATAL_FAILURE(
      {
        auto f = [](const auto& x) { return bad_throws1(x); };

        Eigen::VectorXd x = Eigen::VectorXd::Ones(2);

        stan::test::expect_ad_matvar(f, x);
      },
      "");
}

template <typename T>
auto bad_throws2(const T& x) {
  throw std::exception();
  return x;
}

template <typename T>
auto bad_throws2(const stan::math::var_value<T>& x) {
  return x;
}

TEST(test_unit_math_test_ad_matvar, bad_throws2) {
  EXPECT_FATAL_FAILURE(
      {
        auto f = [](const auto& x) { return bad_throws2(x); };

        Eigen::VectorXd x = Eigen::VectorXd::Ones(2);

        stan::test::expect_ad_matvar(f, x);
      },
      "");
}
