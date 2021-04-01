#ifndef TEST_PROB_UTILITY_HPP
#define TEST_PROB_UTILITY_HPP

#include <stan/math/mix.hpp>

using stan::is_constant_all;
using stan::is_vector;
using stan::scalar_type;
using stan::math::fvar;
using stan::math::var;
using std::vector;

using size_type = stan::math::index_type_t<Eigen::Matrix<double, 1, 1>>;

// ------------------------------------------------------------

struct empty {};

template <typename T>
struct is_empty {
  enum { value = false };
};

template <>
struct is_empty<empty> {
  enum { value = true };
};

// ------------------------------------------------------------

namespace std {
std::ostream& operator<<(std::ostream& os, const vector<double>& param) {
  os << "(";
  for (size_t n = 0; n < param.size(); n++) {
    os << param[n];
    if (n < param.size() - 1)
      os << ", ";
  }
  os << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const vector<var>& param) {
  os << "(";
  for (size_t n = 0; n < param.size(); n++) {
    os << param[n];
    if (n < param.size() - 1)
      os << ", ";
  }
  os << ")";
  return os;
}

}  // namespace std

// ------------------------------------------------------------

template <typename T>
T get_param(const vector<double>& params, const size_t n) {
  T param = 0;
  if (n < params.size())
    param = params[n];
  return param;
}

template <>
empty get_param<empty>(const vector<double>& /*params*/, const size_t /*n*/) {
  return empty();
}

template <>
fvar<double> get_param<fvar<double>>(const vector<double>& params,
                                     const size_t n) {
  fvar<double> param = 0;
  if (n < params.size()) {
    param = params[n];
    param.d_ = 1.0;
  }
  return param;
}
template <>
fvar<var> get_param<fvar<var>>(const vector<double>& params, const size_t n) {
  fvar<var> param = 0;
  if (n < params.size()) {
    param = params[n];
    param.d_ = 1.0;
  }
  return param;
}
template <>
fvar<fvar<double>> get_param<fvar<fvar<double>>>(const vector<double>& params,
                                                 const size_t n) {
  fvar<fvar<double>> param = 0;
  if (n < params.size()) {
    param = params[n];
    param.d_.val_ = 1.0;
  }
  return param;
}
template <>
fvar<fvar<var>> get_param<fvar<fvar<var>>>(const vector<double>& params,
                                           const size_t n) {
  fvar<fvar<var>> param = 0;
  if (n < params.size()) {
    param = params[n];
    param.d_.val_ = 1.0;
  }
  return param;
}

// ------------------------------------------------------------

// default template handles Eigen::Matrix
template <typename T, stan::require_not_var_matrix_t<T>* = nullptr>
T get_params(const vector<vector<double>>& parameters, const size_t p) {
  T param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size())
      param(n) = get_param<stan::scalar_type_t<T>>(parameters[n], p);
  return param;
}

// handle `var_value<T>` where T is an Eigen type
template <typename T, stan::require_var_matrix_t<T>* = nullptr>
T get_params(const vector<vector<double>>& parameters, const size_t p) {
  typename T::value_type param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size())
      param(n) = get_param<double>(parameters[n], p);
  return param;
}

// handle empty
template <>
empty get_params<empty>(const vector<vector<double>>& /*parameters*/,
                        const size_t /*p*/) {
  return empty();
}
// handle scalars
template <>
double get_params<double>(const vector<vector<double>>& parameters,
                          const size_t p) {
  double param(0);
  if (p < parameters[0].size())
    param = parameters[0][p];
  return param;
}
template <>
var get_params<var>(const vector<vector<double>>& parameters, const size_t p) {
  var param(0);
  if (p < parameters[0].size())
    param = parameters[0][p];
  return param;
}
template <>
fvar<double> get_params<fvar<double>>(const vector<vector<double>>& parameters,
                                      const size_t p) {
  fvar<double> param(0);
  if (p < parameters[0].size()) {
    param = parameters[0][p];
    param.d_ = 1.0;
  }
  return param;
}
template <>
fvar<var> get_params<fvar<var>>(const vector<vector<double>>& parameters,
                                const size_t p) {
  fvar<var> param(0);
  if (p < parameters[0].size()) {
    param = parameters[0][p];
    param.d_ = 1.0;
  }
  return param;
}
template <>
fvar<fvar<double>> get_params<fvar<fvar<double>>>(
    const vector<vector<double>>& parameters, const size_t p) {
  fvar<fvar<double>> param(0);
  if (p < parameters[0].size()) {
    param = parameters[0][p];
    param.d_.val_ = 1.0;
  }
  return param;
}
template <>
fvar<fvar<var>> get_params<fvar<fvar<var>>>(
    const vector<vector<double>>& parameters, const size_t p) {
  fvar<fvar<var>> param(0);
  if (p < parameters[0].size()) {
    param = parameters[0][p];
    param.d_.val_ = 1.0;
  }
  return param;
}
template <>
int get_params<int>(const vector<vector<double>>& parameters, const size_t p) {
  int param(0);
  if (p < parameters[0].size())
    param = (int)parameters[0][p];
  return param;
}
// handle vectors
template <>
vector<int> get_params<vector<int>>(const vector<vector<double>>& parameters,
                                    const size_t p) {
  vector<int> param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size())
      param[n] = parameters[n][p];
  return param;
}
template <>
vector<double> get_params<vector<double>>(
    const vector<vector<double>>& parameters, const size_t p) {
  vector<double> param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size())
      param[n] = parameters[n][p];
  return param;
}
template <>
vector<var> get_params<vector<var>>(const vector<vector<double>>& parameters,
                                    const size_t p) {
  vector<var> param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size())
      param[n] = parameters[n][p];
  return param;
}
template <>
vector<fvar<double>> get_params<vector<fvar<double>>>(
    const vector<vector<double>>& parameters, const size_t p) {
  vector<fvar<double>> param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size()) {
      param[n] = parameters[n][p];
      param[n].d_ = 1.0;
    }
  return param;
}
template <>
vector<fvar<var>> get_params<vector<fvar<var>>>(
    const vector<vector<double>>& parameters, const size_t p) {
  vector<fvar<var>> param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size()) {
      param[n] = parameters[n][p];
      param[n].d_ = 1.0;
    }
  return param;
}
template <>
vector<fvar<fvar<double>>> get_params<vector<fvar<fvar<double>>>>(
    const vector<vector<double>>& parameters, const size_t p) {
  vector<fvar<fvar<double>>> param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size()) {
      param[n] = parameters[n][p];
      param[n].d_.val_ = 1.0;
    }
  return param;
}
template <>
vector<fvar<fvar<var>>> get_params<vector<fvar<fvar<var>>>>(
    const vector<vector<double>>& parameters, const size_t p) {
  vector<fvar<fvar<var>>> param(parameters.size());
  for (size_t n = 0; n < parameters.size(); n++)
    if (p < parameters[0].size()) {
      param[n] = parameters[n][p];
      param[n].d_.val_ = 1.0;
    }
  return param;
}

// ------------------------------------------------------------

// default template handles Eigen::Matrix
template <typename T, stan::require_not_var_matrix_t<T>* = nullptr>
T get_params(const vector<vector<double>>& parameters, const size_t /*n*/,
             const size_t p) {
  T param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size())
      param(i) = get_param<stan::scalar_type_t<T>>(parameters[i], p);
  return param;
}

// handle `var_value<T>` where T is an Eigen type
template <typename T, stan::require_var_matrix_t<T>* = nullptr>
T get_params(const vector<vector<double>>& parameters, const size_t /*n*/,
             const size_t p) {
  typename T::value_type param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size())
      param(i) = get_param<double>(parameters[i], p);
  return param;
}

// handle empty
template <>
empty get_params<empty>(const vector<vector<double>>& /*parameters*/,
                        const size_t /*n*/, const size_t /*p*/) {
  return empty();
}
// handle scalars
template <>
double get_params<double>(const vector<vector<double>>& parameters,
                          const size_t n, const size_t p) {
  double param(0);
  if (p < parameters[0].size())
    param = parameters[n][p];
  return param;
}
template <>
var get_params<var>(const vector<vector<double>>& parameters, const size_t n,
                    const size_t p) {
  var param(0);
  if (p < parameters[0].size())
    param = parameters[n][p];
  return param;
}
template <>
fvar<double> get_params<fvar<double>>(const vector<vector<double>>& parameters,
                                      const size_t n, const size_t p) {
  fvar<double> param(0);
  if (p < parameters[0].size()) {
    param = parameters[n][p];
    param.d_ = 1.0;
  }
  return param;
}
template <>
fvar<var> get_params<fvar<var>>(const vector<vector<double>>& parameters,
                                const size_t n, const size_t p) {
  fvar<var> param(0);
  if (p < parameters[0].size()) {
    param = parameters[n][p];
    param.d_ = 1.0;
  }
  return param;
}
template <>
fvar<fvar<double>> get_params<fvar<fvar<double>>>(
    const vector<vector<double>>& parameters, const size_t n, const size_t p) {
  fvar<fvar<double>> param(0);
  if (p < parameters[0].size()) {
    param = parameters[n][p];
    param.d_.val_ = 1.0;
  }
  return param;
}
template <>
fvar<fvar<var>> get_params<fvar<fvar<var>>>(
    const vector<vector<double>>& parameters, const size_t n, const size_t p) {
  fvar<fvar<var>> param(0);
  if (p < parameters[0].size()) {
    param = parameters[n][p];
    param.d_.val_ = 1.0;
  }
  return param;
}
template <>
int get_params<int>(const vector<vector<double>>& parameters, const size_t n,
                    const size_t p) {
  int param(0);
  if (p < parameters[0].size())
    param = (int)parameters[n][p];
  return param;
}
// handle vectors
template <>
vector<int> get_params<vector<int>>(const vector<vector<double>>& parameters,
                                    const size_t /*n*/, const size_t p) {
  vector<int> param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size())
      param[i] = parameters[i][p];
  return param;
}
template <>
vector<double> get_params<vector<double>>(
    const vector<vector<double>>& parameters, const size_t /*n*/,
    const size_t p) {
  vector<double> param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size())
      param[i] = parameters[i][p];
  return param;
}
template <>
vector<var> get_params<vector<var>>(const vector<vector<double>>& parameters,
                                    const size_t /*n*/, const size_t p) {
  vector<var> param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size())
      param[i] = parameters[i][p];
  return param;
}
template <>
vector<fvar<double>> get_params<vector<fvar<double>>>(
    const vector<vector<double>>& parameters, const size_t /*n*/,
    const size_t p) {
  vector<fvar<double>> param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size()) {
      param[i] = parameters[i][p];
      param[i].d_ = 1.0;
    }
  return param;
}
template <>
vector<fvar<var>> get_params<vector<fvar<var>>>(
    const vector<vector<double>>& parameters, const size_t /*n*/,
    const size_t p) {
  vector<fvar<var>> param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size()) {
      param[i] = parameters[i][p];
      param[i].d_ = 1.0;
    }
  return param;
}
template <>
vector<fvar<fvar<double>>> get_params<vector<fvar<fvar<double>>>>(
    const vector<vector<double>>& parameters, const size_t /*n*/,
    const size_t p) {
  vector<fvar<fvar<double>>> param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size()) {
      param[i] = parameters[i][p];
      param[i].d_.val_ = 1.0;
    }
  return param;
}
template <>
vector<fvar<fvar<var>>> get_params<vector<fvar<fvar<var>>>>(
    const vector<vector<double>>& parameters, const size_t /*n*/,
    const size_t p) {
  vector<fvar<fvar<var>>> param(parameters.size());
  for (size_t i = 0; i < parameters.size(); i++)
    if (p < parameters[0].size()) {
      param[i] = parameters[i][p];
      param[i].d_.val_ = 1.0;
    }
  return param;
}

// ------------------------------------------------------------
// default template handles Eigen::Matrix
template <typename T, stan::require_eigen_t<T>* = nullptr>
T get_repeated_params(const vector<double>& parameters, const size_t p,
                      const size_t N_REPEAT) {
  T params(N_REPEAT);
  stan::value_type_t<T> param;

  if (p < parameters.size())
    param = get_param<stan::scalar_type_t<T>>(parameters, p);
  else
    param = 0;

  for (size_t n = 0; n < N_REPEAT; n++) {
    params(n) = param;
  }

  return params;
}

// handle `var_value<T>` where T is an Eigen type
template <typename T, stan::require_var_matrix_t<T>* = nullptr>
T get_repeated_params(const vector<double>& parameters, const size_t p,
                      const size_t N_REPEAT) {
  typename T::value_type params(N_REPEAT);
  double param;

  if (p < parameters.size())
    param = get_param<double>(parameters, p);
  else
    param = 0;

  for (size_t n = 0; n < N_REPEAT; n++) {
    params(n) = param;
  }

  return params;
}

// handle `std::vector`
template <typename T, stan::require_std_vector_t<T>* = nullptr>
T get_repeated_params(const vector<double>& parameters, const size_t p,
                      const size_t N_REPEAT) {
  T params(N_REPEAT);
  stan::value_type_t<T> param;

  if (p < parameters.size())
    param = get_param<stan::scalar_type_t<T>>(parameters, p);
  else
    param = 0;

  for (size_t n = 0; n < N_REPEAT; n++) {
    params[n] = param;
  }

  return params;
}

// handle empty
template <typename T,
          std::enable_if_t<std::is_same<T, empty>::value>* = nullptr>
T get_repeated_params(const vector<double>&, const size_t, const size_t) {
  return T();
}

// handle scalars
template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
T get_repeated_params(const vector<double>& parameters, const size_t p,
                      const size_t /*N_REPEAT*/) {
  if (p < parameters.size())
    return get_param<T>(parameters, p);
  else
    return 0;
}

// ------------------------------------------------------------

template <typename T>
typename scalar_type<T>::type select_var_param(
    const vector<vector<double>>& parameters, const size_t n, const size_t p) {
  typename scalar_type<T>::type param(0);
  if (p < parameters[0].size()) {
    if (!is_constant_all<T>::value)
      param = parameters[n][p];
    else
      param = parameters[0][p];
  }
  return param;
}

template <>
empty select_var_param<empty>(const vector<vector<double>>& /*parameters*/,
                              const size_t /*n*/, const size_t /*p*/) {
  return empty();
}

// ------------------------------------------------------------
template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5>
struct all_scalar {
  enum {
    value = (!is_vector<T0>::value || is_empty<T0>::value)
            && (!is_vector<T1>::value || is_empty<T1>::value)
            && (!is_vector<T2>::value || is_empty<T2>::value)
            && (!is_vector<T3>::value || is_empty<T3>::value)
            && (!is_vector<T4>::value || is_empty<T4>::value)
            && (!is_vector<T5>::value || is_empty<T5>::value)
  };
};

template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5>
struct all_constant {
  enum {
    value = (is_constant_all<T0>::value || is_empty<T0>::value)
            && (is_constant_all<T1>::value || is_empty<T1>::value)
            && (is_constant_all<T2>::value || is_empty<T2>::value)
            && (is_constant_all<T3>::value || is_empty<T3>::value)
            && (is_constant_all<T4>::value || is_empty<T4>::value)
            && (is_constant_all<T5>::value || is_empty<T5>::value)
  };
};

template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5>
struct all_var {
  enum {
    value = (!is_constant_all<T0>::value || is_empty<T0>::value)
            && (!is_constant_all<T1>::value || is_empty<T1>::value)
            && (!is_constant_all<T2>::value || is_empty<T2>::value)
            && (!is_constant_all<T3>::value || is_empty<T3>::value)
            && (!is_constant_all<T4>::value || is_empty<T4>::value)
            && (!is_constant_all<T5>::value || is_empty<T5>::value)
  };
};

template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5>
struct any_vector {
  enum {
    value = is_vector<T0>::value || is_vector<T1>::value || is_vector<T2>::value
            || is_vector<T3>::value || is_vector<T4>::value
            || is_vector<T5>::value
  };
};

// ------------------------------------------------------------
template <typename T, stan::require_not_var_matrix_t<T>* = nullptr>
void add_adjoints(vector<double>& /*x*/, T& /*p*/) {}

template <>
void add_adjoints<var>(vector<double>& x, var& p) {
  x.push_back(p.adj());
}

template <typename T, stan::require_var_matrix_t<T>* = nullptr>
void add_adjoints(vector<double>& x, T& p) {
  for (size_type n = 0; n < p.size(); n++) {
    x.push_back(p.adj().coeff(n));
  }
}

template <>
void add_adjoints<vector<var>>(vector<double>& x, vector<var>& p) {
  for (size_type n = 0; n < p.size(); n++)
    x.push_back(p[n].adj());
}

template <>
void add_adjoints<Eigen::Matrix<var, 1, Eigen::Dynamic>>(
    vector<double>& x, Eigen::Matrix<var, 1, Eigen::Dynamic>& p) {
  for (size_type n = 0; n < p.size(); n++)
    x.push_back(p(n).adj());
}

template <>
void add_adjoints<Eigen::Matrix<var, Eigen::Dynamic, 1>>(
    vector<double>& x, Eigen::Matrix<var, Eigen::Dynamic, 1>& p) {
  for (size_type n = 0; n < p.size(); n++)
    x.push_back(p(n).adj());
}

template <>
void add_adjoints<fvar<var>>(vector<double>& x, fvar<var>& p) {
  x.push_back(p.val_.adj());
}

template <>
void add_adjoints<vector<fvar<var>>>(vector<double>& x, vector<fvar<var>>& p) {
  for (size_t n = 0; n < p.size(); n++)
    x.push_back(p[n].val_.adj());
}

template <>
void add_adjoints<Eigen::Matrix<fvar<var>, 1, Eigen::Dynamic>>(
    vector<double>& x, Eigen::Matrix<fvar<var>, 1, Eigen::Dynamic>& p) {
  for (size_type n = 0; n < p.size(); n++)
    x.push_back(p(n).val_.adj());
}

template <>
void add_adjoints<Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1>>(
    vector<double>& x, Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1>& p) {
  for (size_type n = 0; n < p.size(); n++)
    x.push_back(p(n).val_.adj());
}

template <>
void add_adjoints<fvar<fvar<var>>>(vector<double>& x, fvar<fvar<var>>& p) {
  x.push_back(p.val_.val_.adj());
}

template <>
void add_adjoints<vector<fvar<fvar<var>>>>(vector<double>& x,
                                           vector<fvar<fvar<var>>>& p) {
  for (size_t n = 0; n < p.size(); n++)
    x.push_back(p[n].val_.val_.adj());
}

template <>
void add_adjoints<Eigen::Matrix<fvar<fvar<var>>, 1, Eigen::Dynamic>>(
    vector<double>& x, Eigen::Matrix<fvar<fvar<var>>, 1, Eigen::Dynamic>& p) {
  for (size_type n = 0; n < p.size(); n++)
    x.push_back(p(n).val_.val_.adj());
}

template <>
void add_adjoints<Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1>>(
    vector<double>& x, Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1>& p) {
  for (size_type n = 0; n < p.size(); n++)
    x.push_back(p(n).val_.val_.adj());
}

template <typename T0, typename T1, typename T2, typename T3, typename T4,
          typename T5>
void add_adjoints(vector<double>& x, T0& p0, T1& p1, T2& p2, T3& p3, T4& p4,
                  T5& p5) {
  if (!is_constant_all<T0>::value)
    add_adjoints(x, p0);
  if (!is_constant_all<T1>::value)
    add_adjoints(x, p1);
  if (!is_constant_all<T2>::value)
    add_adjoints(x, p2);
  if (!is_constant_all<T3>::value)
    add_adjoints(x, p3);
  if (!is_constant_all<T4>::value)
    add_adjoints(x, p4);
  if (!is_constant_all<T5>::value)
    add_adjoints(x, p5);
}

#endif
