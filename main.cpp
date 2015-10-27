#include <iostream>

#include <stan/math/rev/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system_cvode.hpp>
#include <stan/math/prim/arr/functor/integrate_ode_cvode.hpp>
#include <stan/math/rev/mat/functor/gradient.hpp>

// clang++ -O3 -I . -isystem lib/eigen_3.2.4 -isystem lib/boost_1.58.0 -Ilib/cvode_2.8.2/include -L/usr/local/lib/ -lsundials_cvode -lsundials_nvecserial -o main main.cpp

class sho_functor {
public:
  template <typename T0, typename T1, typename T2>
  inline
  std::vector<typename stan::return_type<T1, T2>::type>
  operator()(const T0& t_in,                 // time
             const std::vector<T1>& y_in,    // state
             const std::vector<T2>& theta,   // parameters
             const std::vector<double>& x,   // double data
             const std::vector<int>& x_int,  // integer data
             std::ostream* msgs) const {
    if (y_in.size() != 2)
      throw std::domain_error("Functor called with inconsistent state");
    
    std::vector<typename stan::return_type<T1, T2>::type> f;
    f.push_back(y_in.at(1));
    f.push_back(- theta.at(0) * theta.at(0) * y_in.at(0));
    
    return f;
  }
};

class test_functor_double_var_1 {
public:
  template <typename T>
  inline
  T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;
    
    std::vector<T> theta;
    theta.push_back(x(0));
    
    std::vector<double> y0;
    y0.push_back(1.25);
    y0.push_back(0.0);
    
    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);
    
    std::vector<double> data;
    std::vector<int> data_int;
    
    std::vector<std::vector<T> > ys
    = stan::math::integrate_ode_cvode(sho, y0, t0, ts, theta, data, data_int);
    
    return ys[0][0];
  }
};

int main() {

  sho_functor sho;
  
  std::vector<stan::math::var> theta;
  theta.push_back(0.5);
  
  std::vector<double> y0;
  y0.push_back(1.25);
  y0.push_back(0.0);
  
  double t0 = 0.0;
  std::vector<double> ts;
  ts.push_back(5.0);
  
  std::vector<double> data;
  std::vector<int> data_int;
  
  stan::math::coupled_ode_system_cvode<sho_functor, double, stan::math::var>
    coupled_system(sho, y0, t0, theta, data, data_int, 1e-8, 1e-10, 1e6, 0);
  
  std::vector<std::vector<double> > y_coupled(ts.size());
  for (int n = 0; n < ts.size(); ++n)
    y_coupled[n].resize(y0.size());
  
  coupled_system.integrate_times(ts, y_coupled);
  
  /*
  double y[2];
  
  double* J[4];
  for (int n = 0; n < 4; ++n)
    J[n] = new double[2];
  
  long int s_mu = 0;
  double t = 5;
  
  coupled_system.banded_jacobian(y, J, 2, t);

  for (int n = 0; n < 4; ++n)
    delete[] J[n];
  */
  /*
  double omega = 0.5;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x_var(1);
  x_var(0) = omega;
  test_functor_double_var_1 func1;
  func1(x_var);
  */

  return 0;
}
