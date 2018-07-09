#ifdef STAN_OPENSEES

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <algorithm>

#include <iostream>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <StandardStream.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <SectionRepres.h>
#include <TimeSeries.h>
#include <CrdTransf.h>
#include <BeamIntegration.h>
#include <NodalLoad.h>
#include <AnalysisModel.h>
#include <PlainHandler.h>
#include <RCM.h>
#include <AMDNumberer.h>
#include <LimitCurve.h>
#include <DamageModel.h>
#include <FrictionModel.h>
#include <HystereticBackbone.h>
#include <YieldSurface_BC.h>
#include <CyclicModel.h>
#include <FileStream.h>
#include <CTestNormUnbalance.h>
#include <NewtonRaphson.h>
#include <TransformationConstraintHandler.h>
#include <Newmark.h>
#include <ProfileSPDLinSolver.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ProfileSPDLinSOE.h>
#include <SymBandEigenSolver.h>
#include <SymBandEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <ArpackSOE.h>
#include <LoadControl.h>
#include <CTestPFEM.h>
#include <PFEMIntegrator.h>
#include <TransientIntegrator.h>
#include <PFEMSolver.h>
#include <PFEMLinSOE.h>
#include <Accelerator.h>
#include <KrylovAccelerator.h>
#include <AcceleratedNewton.h>
#include <RaphsonAccelerator.h>
#include <SecantAccelerator2.h>
#include <PeriodicAccelerator.h>
#include <LineSearch.h>
#include <InitialInterpolatedLineSearch.h>
#include <BisectionLineSearch.h>
#include <SecantLineSearch.h>
#include <RegulaFalsiLineSearch.h>
#include <NewtonLineSearch.h>
#include <FileDatastore.h>
#include <UniaxialMaterial.h>
#include <Domain.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <EquiSolnAlgo.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBrokerAllClasses.h>
#include <PFEMAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#ifdef _RELIABILITY
#include <ReliabilityStaticAnalysis.h>
#include <ReliabilityDirectIntegrationAnalysis.h>
#endif
#include <Timer.h>
#include <SimulationInformation.h>
#include <SP_Constraint.h>
#include <LinearCrdTransf2d.h>
#include <ElasticBeam2d.h>
#include <LinearSeries.h>
#include <LoadPattern.h>
#include <DOF_Numberer.h>

using namespace std;

class OpenSEESModel {
  
public:
  // inline std::vector<std::vector<double> >
  // operator()(const std::vector<double>& theta,
  //            const bool need_sens,
  //            const std::vector<double>& x_r,
  //            const std::vector<int>& x_i,
  //            std::ostream* msgs = nullptr) const {
  //   std::vector<std::vector<double> > res;
  //   if(need_sens)
  //     res = solve_with_sensitivity(theta);
  //   else
  //     res = solve(theta);

  //   return res;
  // }

  // // sensitivity calc: finite diff
  // inline std::vector<std::vector<double> >
  // solve_with_sensitivity(const std::vector<double>& theta) const {
  //   const double h = 1.e-3;
  //   std::vector<double> theta2{theta[0] + h};
  //   std::vector<std::vector<double> > sol1 = solve(theta);
  //   std::vector<std::vector<double> > sol2 = solve(theta2);
  //   const double sens = (sol2[0][0] - sol1[0][0])/h;
  //   return {{sol1[0][0], sens}};
  // }

  // // QoI-only, as MFEM has no sensitivity module
  // inline std::vector<std::vector<double> >
  // solve(const std::vector<double>& theta) const {
    
  // }
};


TEST(pde_solvers, opensees_sensitivity) {
  using stan::math::forward_pde;

  // OpenSEESModel pde;
  // std::vector<double> x_r;
  // std::vector<int> x_i;
  
  // std::vector<double> theta{1.0};
  // const double err_exact = 0.015791582;
  // std::vector<double> qoi = forward_pde(pde, theta, x_r, x_i);
  // ASSERT_FLOAT_EQ(qoi[0], err_exact);

  // std::vector<stan::math::var> theta_v{1.2};
  // std::vector<stan::math::var> qoi_v = forward_pde(pde, theta_v, x_r, x_i);
  // std::vector<double> g;
  // stan::math::set_zero_all_adjoints();
  // qoi_v[0].grad(theta_v, g);
  // const double sens_exact = 2.1105776;
  // ASSERT_FLOAT_EQ(g[0], sens_exact);
}

#endif
