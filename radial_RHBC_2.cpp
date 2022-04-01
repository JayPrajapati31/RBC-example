#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <rbc_ProfileCuttingManager.h>
#include <rbc_ProfileCuttingManager_impl.h>
#include <rbc_CubicBezierCurve.h>
#include <rbc_CubicSplineCurve.h>
#include <rbc_ParametricCurveCombination.h>
#include <math.h>
#include "customFunc.h"

void weightFuncFixNode(const rbc::WireConfig config, const double* weightRange, const rbc::ParametricCurve& init, int& node, double& tparam);
void weightFuncFixNode(const int nNodes,const double blade_len,const double* bounds, const double* weightRange, const rbc::ParametricCurve& init,  int& node, double& tparam);

int main(int argc, char** argv)
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Create a uniform 1D mesh on the interval [0,1]
  const double kval = 50.;
  const int nNodes = 101;
  const int nElm = nNodes-1;
  double blade_length = 0.5;
  const double block_size = 0.25;
  double redundentLen = 0.2525;
  bool SmallBlade_FineCut = true;

  // Target shape
  rbc::NLopt_Params nlopt_params;
  nlopt_params.set_ftol_abs = true; nlopt_params.set_ftol_rel = true;
  nlopt_params.set_xtol_abs = true; nlopt_params.set_xtol_rel = true;
  const double tol = 1.e-10;
  nlopt_params.ftol_abs = tol;
  nlopt_params.ftol_rel = tol;
  nlopt_params.xtol_abs = tol;
  nlopt_params.xtol_rel = tol;

  rbc::Ipopt_Params ipopt_params;
  ipopt_params.derivative_test = "first-order";
  ipopt_params.obj_scaling_factor = 0.1;
 
  //input target shape
  std::vector<double> xy, tvec;
  double dt = 1.0/100.;
   for(int i=0;i<101;  ++i)
    {
      tvec.push_back(i*dt);
      double x = tvec[i]*block_size;
      double y = 2.*x*(x-0.25) +0.5;
      xy.push_back(x);
      xy.push_back(y);
    }
   assert(tvec.size() == xy.size()/2);
   int nknots = tvec.size(),nSample =20;
   std::vector<double>  init_xy(2*nknots), inter_xy(2*nknots), init_x(nknots), init_y(nknots), target_y(nknots);
   //  readTarget("/home/jay/CodeRepos/rbc/examples/axisymm/Radial RHBC/matlab files/target.txt", xy,  tvec);

  for(int i=0; i<nknots; ++i)
    {
      target_y[i] = xy[2*i+1];
      init_xy[2*i] = xy[2*i];
      init_xy[2*i+1] = redundentLen;  // (x,y) -> straight curve
      init_x[i] = init_xy[2*i];
      init_y[i] = init_xy[2*i+1];
    }

  rbc::CubicSplineCurve target_curve(tvec,nSample,nlopt_params), init_curve(tvec,nSample,nlopt_params), intermediate_curve(tvec,nSample,nlopt_params), buffer_curve(tvec,nSample, nlopt_params);
  target_curve.SetControlPoints(xy);
  init_curve.SetControlPoints(init_xy);

  double y_lBound, y_uBound,x_lBound, x_uBound, SmallBlade_FineCut_y0;
  y_lBound = *min_element(init_y.begin(), init_y.end()); // min ycoordinate --> y coordinate od straight initial shape
  y_uBound = *max_element(target_y.begin(), target_y.end()); // max y coordinate --> max y coordinate of target shape
  x_lBound = *min_element(init_x.begin(), init_x.end()); 
  x_uBound = *max_element(init_x.begin(), init_x.end());
  SmallBlade_FineCut_y0 = *min_element(target_y.begin(), target_y.end()) - 0.15 ;
  
  // linear combination of profiles
  rbc::ParametricCurveCombination curve_combination(target_curve,init_curve,100,nlopt_params);

  // Profile cutting manager
  rbc::ProfileCuttingManager manager(target_curve, std::vector<double>{0.5*(x_lBound+x_uBound) - 0.5*block_size,redundentLen,
									 0.5*(x_lBound+x_uBound) - 0.5*block_size, y_uBound + 0.1,
									 0.5*(x_lBound+x_uBound) + 0.5*block_size, y_uBound+0.1,
									 0.5*(x_lBound+x_uBound) + 0.5*block_size, redundentLen});
  manager.CreateBlade(blade_length, nElm);
  
  // Constants
  const double sdTOL = 1.e-6;  // determines which points on the target has been cut
  const double dtTOL = 0.01; // size of smallest cut = size of 4 element
    
  // Solve parameters
  rbc::SolveParams solve_params({.EPS = 1.e-6, .resScale = 1., .dofScale = 1., .nMaxIters = 25, .verbose = false});
  
  // define optimization variables
  std::array<bool,6> optim_vars{true, true, true, true, false, true}; // H, V, M, BC, r=y
  
  // bounds for optimization parameters
  std::vector<double> lbounds(5), ubounds(5),ip_lbounds(5), ip_ubounds(5) ;
  lbounds[0] = -HUGE_VAL;//100.2;
  lbounds[1] = -HUGE_VAL;//100.2;
  lbounds[2] = -HUGE_VAL;//2.*(M_PI/blade_length);
  lbounds[3] = -(M_PI/180.)*90.;
  lbounds[4] =  y_lBound-0.1;
  ip_lbounds[4] = lbounds[4];
  
  ubounds[0] = HUGE_VAL;//100.2;
  ubounds[1] = HUGE_VAL;//100.2;
  ubounds[2] = HUGE_VAL;//2.0*(M_PI/blade_length);
  ubounds[3] = (M_PI/180.)*90.;
  ubounds[4] = y_uBound + 0.1;
  ip_ubounds[4] = ubounds[4];
  
  // path to store results
  std::string op_path = "/home/jay/CodeRepos/rbc/examples/axisymm/Radial RHBC/Data/"; 
  
  // alpha and deltaAlpha
  int m = 11;  // # of intermeadiate profiles
  int opt_count = 0;
  double dAlpha = (1.0)/static_cast<double>(m-1);
  double t_fix = 0.5;
  const double advance_vec[] = {0., 1.0};
  std::vector<std::vector<double>> endData;
  std::vector<double> temp_end_data(6,55);

  // Wire parameters
  rbc::WireParams in_wire_params, out_wire_params, temp_wire_params, nlopt_wire_params;
  in_wire_params.H = 0.0;
  in_wire_params.V = 0.0;
  in_wire_params.M = 0.0;
  in_wire_params.BC = 0.0;
  in_wire_params.node = nNodes/2;
  
  // Wire cofnfiguration
  rbc::WireConfig in_wire_config, out_wire_config, temp_wire_config, nlopt_wire_config;
  in_wire_config.theta = std::vector<double>(nNodes, 0.);
  in_wire_config.xy = std::vector<double>(2*nNodes,0.);

  for(int j=0; j<nNodes; ++j)
    {
      in_wire_config.xy[2*j] = double(j)*static_cast<double>(blade_length/nElm);
      in_wire_config.xy[2*j+1] = redundentLen;
    }
  
  // looping over intermediate shapes
  double X[2]; int feasible = 0;
  double box_min_y, elastica_max_y;
  for(int profile=0; profile<m; ++profile)
    {      
      std::cout<< "\n*************************\n"<< "\n\nOptimizing profile # "<<profile<< "\n\n"<< std::endl;     
      curve_combination.SetCombinationFactor(profile*dAlpha); // setting alpha
      for(int j=0; j<nknots; ++j)	
	{
	  curve_combination.Evaluate(tvec[j],X);
	  inter_xy[2*j] = X[0];
	  inter_xy[2*j+1] = X[1];
	}
      intermediate_curve.SetControlPoints(inter_xy);
      std::cout << "\n\n t_fix : " << t_fix << std::endl;
      
      //Setting buffer profile as second last profile instead of bounding box for cut region
      if(profile == m-2)
	buffer_curve.SetControlPoints(inter_xy);
      
      PrintFile(op_path+"profile-"+std::to_string(profile)+".dat", intermediate_curve);     
      
      // midpoint 
      intermediate_curve.Evaluate(t_fix, in_wire_params.xy, nullptr);
      if(profile == 0)
	in_wire_params.xy[0] -= 0.01;
      
      // Advance the cut
      manager.AdvanceCut(optim_vars, in_wire_params, lbounds, ubounds, in_wire_config, intermediate_curve,
			 solve_params, ipopt_params, sdTOL, dtTOL, advance_vec, out_wire_params, out_wire_config, 0.0, 0.0);
      
      std::cout << "\nis wire feasible : " <<  manager.IsWireFeasible(out_wire_config,sdTOL) << std::endl;
     
      // Output  = input for next cutting step
      in_wire_config.theta = out_wire_config.theta;
      in_wire_config.xy = out_wire_config.xy;
      in_wire_params = out_wire_params;

      PrintOpt(op_path+"opt-"+std::to_string(opt_count)+".dat", manager);
      CheckConstraint(target_curve,manager.GetBlade());
      
      //printig remaining area points
      PrintFile(op_path+"rem-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
      box_min_y = PrintBox(op_path+"box-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
      double temp_xy[2];
      int temp_node;
      temp_node =  ElasticaMax(manager.GetBlade(), temp_xy);
      elastica_max_y = temp_xy[1];
      PrintIntervals(manager.GetCutProfileMonitor());
      
      temp_end_data[0] = (out_wire_config.xy)[0];
      temp_end_data[1] = (out_wire_config.xy)[1];
      temp_end_data[2] = (out_wire_config.theta)[0];
      temp_end_data[3] = (out_wire_config.xy[2*(nNodes-1)]);
      temp_end_data[4] = (out_wire_config.xy[2*(nNodes-1)+1]);   
      temp_end_data[5] = (out_wire_config.theta[nNodes-1]);      
      endData.push_back(temp_end_data);
      ++opt_count;
      
      if(profile == m-1)
	{
	  temp_wire_params.H = out_wire_params.H;
	  temp_wire_params.V = out_wire_params.V;
	  temp_wire_params.M = out_wire_params.M;
	  temp_wire_params.BC = out_wire_params.BC;
	  temp_wire_params.xy[0] = out_wire_params.xy[0];
	  temp_wire_params.xy[1] = out_wire_params.xy[1];
	  temp_wire_params.node = out_wire_params.node;
	  temp_wire_config = out_wire_config;
	}
    }//rough cut loop in profiles ends here
  
  std::cout << "\n\n finer cuts " << std::endl;
  optim_vars[5] = true;
  rbc::WeightFunction wfunc(kval);
  int nRoughCuts_longBlade = 2;
  double uncutRange[2], xBound[2] = {x_lBound, x_uBound};

  for(int nIter =0; nIter<nRoughCuts_longBlade; ++nIter)
    {
      std::vector<double> uncut = uncutVec(manager.GetCutProfileMonitor());
      for(int k=0; k<(int)(uncut.size()/2); ++k)
	{
	  std::cout << "\n\n fine cuts # \t" << k << std::endl;
	  wfunc.SetIntervals(std::vector<std::pair<double,double>>{std::make_pair( uncut[2*k] - 0.05, uncut[2*k+1]+0.05)});
	  uncutRange[0] = uncut[2*k];
	  uncutRange[1] = uncut[2*k+1];
	  weightFuncFixNode(nNodes,blade_length, xBound ,uncutRange,init_curve, temp_wire_params.node, t_fix);
	  target_curve.Evaluate(t_fix, temp_wire_params.xy);

	  manager.AdvanceCut(optim_vars, temp_wire_params, lbounds, ubounds, temp_wire_config, target_curve,wfunc,			     
			     solve_params, nlopt_params, sdTOL, dtTOL, advance_vec, out_wire_params, out_wire_config, 0.0, 0.0);
	    
	  double temp_xy[2];
	  int temp_node;
	  temp_node =  ElasticaMax(manager.GetBlade(), temp_xy);
	  elastica_max_y = temp_xy[1];
	  PrintOpt(op_path+"opt-"+std::to_string(opt_count)+".dat", manager);
	  PrintFile(op_path+"rem-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
	  box_min_y = PrintBox(op_path+"box-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
	  PrintIntervals(manager.GetCutProfileMonitor());

	  temp_end_data[0] = (out_wire_config.xy)[0];
	  temp_end_data[1] = (out_wire_config.xy)[1];
	  temp_end_data[2] = (out_wire_config.theta)[0];
	  temp_end_data[3] = (out_wire_config.xy[2*nNodes-2]);
	  temp_end_data[4] = (out_wire_config.xy[2*nNodes-1]);   
	  temp_end_data[5] = (out_wire_config.theta[nNodes-1]);      
	  endData.push_back(temp_end_data);
	  ++opt_count;
	}
      uncut.clear();
    }
}// main function ends

// fix nodes on the target profile
void weightFuncFixNode(const int nNodes,const double blade_len,const double* bounds, const double* weightRange, const rbc::ParametricCurve& init,  int& node, double& tparam)
{
  const double x0 = 0.5*(bounds[0] + bounds[1]) - 0.5*blade_len;
  const double xl = 0.5*(bounds[0] + bounds[1]) + 0.5*blade_len;
  const double dx = (xl-x0)/(nNodes-1);
  std::vector<double> x_coord(nNodes,0);
  
  for(int i=0; i<nNodes; ++i)
    {
      x_coord[i] = x0 + i*dx;
    }

  double t_mid = 0.5*(weightRange[0]+weightRange[1]);
  double x_mid = t_mid*(bounds[1]-bounds[0]);
  int node0, node1;

  for(int i=0; i<nNodes-1; ++i)
    {
      if(x_mid >= x_coord[i] && x_mid <= x_coord[i+1])
	{
	  node = i;
	}
    }

  tparam = (x_coord[node] - bounds[0])/(bounds[1]-bounds[0]);
  if(tparam <0)
    tparam = 0;
  if(tparam >1)
    tparam = 1;
  
  std::cout << "\n\nnode \t " << node << "\t" << x_coord[node] << "param \t" << tparam << std::endl;
}
