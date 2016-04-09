#pragma once

#include <stan/math/rev/arr/functor/ode_model.hpp>

#include "hornberg.hpp"

// define an analytic jacobian for the hornberg for Stan's CVODES
// solver by partial template specialisation


namespace stan {
  namespace math {

    // todo: idealy, we derive the specialised ode_model from the
    // general one which gives possibly syntax problem which is why we
    // chose to use a typedef, maybe we need another helper... not
    // sure.
    template<>
    struct ode_model<hornberg_ode_fun> {
      typedef hornberg_ode_fun F;
      const F& f_;
      const std::vector<double>& theta_;
      const std::vector<double>& x_;
      const std::vector<int>& x_int_;
      std::ostream* msgs_;

      ode_model(const F& f,
		const std::vector<double>& theta,
		const std::vector<double>& x,
		const std::vector<int>& x_int,
		std::ostream* msgs)
	: f_(f),
	  theta_(theta),
	  x_(x),
	  x_int_(x_int),
	  msgs_(msgs)
      {}

      void operator()(const std::vector<double>& y,
		      std::vector<double>& dy_dt,
		      const double t) const {
	dy_dt = f_(t, y, theta_, x_, x_int_, msgs_);
      }

      template <typename Derived1, typename Derived2>
      void
      jacobian_S(const double t,
		 const std::vector<double>& y,
		 Eigen::MatrixBase<Derived1>& fy,
		 Eigen::MatrixBase<Derived2>& Jy
		 ) const {
	using Eigen::VectorXd;
	using std::vector;
	using std::pow;

	const vector<double> f = f_(t, y, theta_, x_, x_int_, msgs_);
	fy = VectorXd::Map(&f[0], f.size());

	const double R = y[0];
	const double Rin = y[1];
	const double x1 = y[2];
	const double x1p = y[3];
	const double x2 = y[4];
	const double x2p = y[5];
	const double x3 = y[6];
	const double x3p = y[7];

	const vector<double>& parms = theta_;

	Jy.setZero();
	Eigen::Map<Eigen::VectorXd> J = Eigen::Map<Eigen::VectorXd>(&Jy(0,0), Jy.cols()*Jy.rows());

	J[0] = -(Vm1/(Km1+R)-(Vm1*R)/pow((Km1+R),2));
	J[1] = Vm1/(Km1+R)-(Vm1*R)/pow((Km1+R),2);
	J[2] = -(k3*x1/(Km3+x1));
	J[3] = k3*x1/(Km3+x1);
	J[8] = Vm2/(Km2+Rin)-(Vm2*Rin)/pow((Km2+Rin),2);
	J[9] = -(Vm2/(Km2+Rin)-(Vm2*Rin)/pow((Km2+Rin),2));
	J[18] = -(k3*R/(Km3+x1)-(k3*R*x1)/pow((Km3+x1),2));
	J[19] = k3*R/(Km3+x1)-(k3*R*x1)/pow((Km3+x1),2);
	J[26] = Vm4/(Km4+x1p)-(Vm4*x1p)/pow((Km4+x1p),2);
	J[27] = -(Vm4/(Km4+x1p)-(Vm4*x1p)/pow((Km4+x1p),2));
	J[28] = -(k5*x2/(Km5+x2));
	J[29] = k5*x2/(Km5+x2);
	J[36] = -(k5*x1p/(Km5+x2)-(k5*x1p*x2)/pow((Km5+x2),2));
	J[37] = k5*x1p/(Km5+x2)-(k5*x1p*x2)/pow((Km5+x2),2);
	J[44] = Vm6/(Km6+x2p)-(Vm6*x2p)/pow((Km6+x2p),2);
	J[45] = -(Vm6/(Km6+x2p)-(Vm6*x2p)/pow((Km6+x2p),2));
	J[46] = -(k7*x3/(Km7+x3));
	J[47] = k7*x3/(Km7+x3);
	J[54] = -(k7*x2p/(Km7+x3)-(k7*x2p*x3)/pow((Km7+x3),2));
	J[55] = k7*x2p/(Km7+x3)-(k7*x2p*x3)/pow((Km7+x3),2);
	J[62] = Vm8/(Km8*(1+Inh/Ki8+x3p/Km8))-(Vm8*x3p)*(Km8*(1/Km8))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2);
	J[63] = -(Vm8/(Km8*(1+Inh/Ki8+x3p/Km8))-(Vm8*x3p)*(Km8*(1/Km8))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2));
      }

      template <typename Derived1, typename Derived2, typename Derived3>
      void
      jacobian_SP(const double t,
		  const std::vector<double>& y,
		  Eigen::MatrixBase<Derived1>& fy,
		  Eigen::MatrixBase<Derived2>& Jy,
		  Eigen::MatrixBase<Derived3>& Jtheta
		  ) const {
	using Eigen::VectorXd;
	using std::vector;
	using std::pow;

	jacobian_S(t, y, fy, Jy);

	const double R = y[0];
	const double Rin = y[1];
	const double x1 = y[2];
	const double x1p = y[3];
	const double x2 = y[4];
	const double x2p = y[5];
	const double x3 = y[6];
	const double x3p = y[7];

	const vector<double>& parms = theta_;

	Jtheta.setZero();

	Eigen::Map<Eigen::VectorXd> J = Eigen::Map<Eigen::VectorXd>(&Jtheta(0,0), Jtheta.cols()*Jtheta.rows());

	J[0] = -(R/(Km1+R));
	J[1] = R/(Km1+R);
	J[8] = (Vm1*R)/pow((Km1+R),2);
	J[9] = -((Vm1*R)/pow((Km1+R),2));
	J[16] = Rin/(Km2+Rin);
	J[17] = -(Rin/(Km2+Rin));
	J[24] = -((Vm2*Rin)/pow((Km2+Rin),2));
	J[25] = (Vm2*Rin)/pow((Km2+Rin),2);
	J[34] = -(R*x1/(Km3+x1));
	J[35] = R*x1/(Km3+x1);
	J[42] = (k3*R*x1)/pow((Km3+x1),2);
	J[43] = -((k3*R*x1)/pow((Km3+x1),2));
	J[50] = x1p/(Km4+x1p);
	J[51] = -(x1p/(Km4+x1p));
	J[58] = -((Vm4*x1p)/pow((Km4+x1p),2));
	J[59] = (Vm4*x1p)/pow((Km4+x1p),2);
	J[68] = -(x1p*x2/(Km5+x2));
	J[69] = x1p*x2/(Km5+x2);
	J[76] = (k5*x1p*x2)/pow((Km5+x2),2);
	J[77] = -((k5*x1p*x2)/pow((Km5+x2),2));
	J[84] = x2p/(Km6+x2p);
	J[85] = -(x2p/(Km6+x2p));
	J[92] = -((Vm6*x2p)/pow((Km6+x2p),2));
	J[93] = (Vm6*x2p)/pow((Km6+x2p),2);
	J[102] = -(x2p*x3/(Km7+x3));
	J[103] = x2p*x3/(Km7+x3);
	J[110] = (k7*x2p*x3)/pow((Km7+x3),2);
	J[111] = -((k7*x2p*x3)/pow((Km7+x3),2));
	J[118] = x3p/(Km8*(1+Inh/Ki8+x3p/Km8));
	J[119] = -(x3p/(Km8*(1+Inh/Ki8+x3p/Km8)));
	J[126] = -((Vm8*x3p)*((1+Inh/Ki8+x3p/Km8)-Km8*(x3p/pow(Km8,2)))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2));
	J[127] = (Vm8*x3p)*((1+Inh/Ki8+x3p/Km8)-Km8*(x3p/pow(Km8,2)))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2);
	J[134] = -((Vm8*x3p)*(Km8*(1/Ki8))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2));
	J[135] = (Vm8*x3p)*(Km8*(1/Ki8))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2);
	J[142] = (Vm8*x3p)*(Km8*(Inh/pow(Ki8,2)))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2);
	J[143] = -((Vm8*x3p)*(Km8*(Inh/pow(Ki8,2)))/pow((Km8*(1+Inh/Ki8+x3p/Km8)),2));
      }

    };

  }
}
