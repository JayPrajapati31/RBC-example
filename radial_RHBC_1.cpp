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

void printTar(std::string fname, const rbc::ParametricCurve& curve);
void weightFuncFixNode(const int nNodes,const double blade_len,const double* bounds, const double* weightRange, const rbc::ParametricCurve& init,  int& node, double& tparam);

int main(int argc, char** argv)
{
  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Create a uniform 1D mesh on the interval [0,1]
  const double kval = 100.;
  const int nNodes = 101;
  const int nElm = nNodes-1;
  double blade_length = 0.5;
  double redundentLen = 0.2525;
  bool SmallBlade_FineCut = true;

  // Target shape
  rbc::NLopt_Params nlopt_params;
  nlopt_params.set_ftol_abs = false ;nlopt_params.set_ftol_rel = true;
  nlopt_params.set_xtol_abs = false; nlopt_params.set_xtol_rel = true;
  const double ftol = 1.e-6;
  const double xtol = 1.e-12;
  nlopt_params.ftol_abs = ftol;
  nlopt_params.ftol_rel = ftol;
  nlopt_params.xtol_abs = xtol;
  nlopt_params.xtol_rel = xtol;

  rbc::Ipopt_Params ipopt_params; 

  std::vector<double> tvec, xy;
  readTarget("../target.txt", xy, tvec);
  assert(tvec.size() == xy.size()/2);
  int nknots = tvec.size();
  int nSample = 20;
  std::vector<double> init_xy(2*nknots), inter_xy(2*nknots), init_y(nknots), target_y(nknots);
  double y_lBound, y_uBound, SmallBlade_FineCut_y0;
  double block_size = 0.25;//0.25;

  for(int i=0; i<nknots; ++i)
    {    
      target_y[i] = xy[2*i+1];      
      init_xy[2*i] = block_size*tvec[i];
      init_xy[2*i+1] = redundentLen;  // (x,y) -> straight curve      
      init_y[i] = init_xy[2*i+1];
    }
  
  rbc::CubicSplineCurve target_curve(tvec,nSample,nlopt_params), init_curve(tvec,nSample,nlopt_params), intermediate_curve(tvec,nSample,nlopt_params), buffer_curve(tvec,nSample, nlopt_params);
  target_curve.SetControlPoints(xy);
  init_curve.SetControlPoints(init_xy);

  y_lBound = *min_element(init_y.begin(), init_y.end());
  y_uBound = *max_element(target_y.begin(), target_y.end());
  SmallBlade_FineCut_y0 = *min_element(target_y.begin(), target_y.end()) - 0.15 ;
  
  // linear combination of profiles
  rbc::ParametricCurveCombination curve_combination(target_curve,init_curve,100,nlopt_params);

  // Profile cutting manager
  rbc::ProfileCuttingManager manager(target_curve, std::vector<double>{0.,redundentLen,
									 0., y_uBound + 0.1,
									 block_size, y_uBound+0.1,
									 block_size, redundentLen});
  manager.CreateBlade(blade_length, nElm);  
    
  // Constants
  const double sdTOL = 0.001;  // determines which points on the target has been cut
  const double dtTOL = 0.005; // size of smallest cut = size of 4 element
    
  // Solve parameters
  rbc::SolveParams solve_params({.EPS = 1.e-6, .resScale = 1., .dofScale = 1., .nMaxIters = 50, .verbose = false});
  
  // define optimization variables
  std::array<bool,6> optim_vars{true, true, true, true, false, true}; // H, V, M, BC, r=y
  
  // bounds for optimization parameters
  std::vector<double> lbounds(5), ubounds(5),ip_lbounds(5), ip_ubounds(5) ;
  lbounds[0] = -HUGE_VAL;//100.2;
  lbounds[1] = -HUGE_VAL;//100.2;
  lbounds[2] = -HUGE_VAL;//2.0*(M_PI/blade_length);
  lbounds[3] = -(M_PI/180.)*180.;
  lbounds[4] =  y_lBound-0.2;

  ubounds[0] = HUGE_VAL;//100.2;
  ubounds[1] = HUGE_VAL;//100.2;
  ubounds[2] = HUGE_VAL;//2.*(M_PI/blade_length);
  ubounds[3] = (M_PI/180.)*180.;
  ubounds[4] = y_uBound + 0.2;
  
  // path to store results
  std::string op_path = "/home/jay/CodeRepos/rbc/examples/axisymm/Radial RHBC/Data/"; 
  
  // alpha and deltaAlpha
  int m = 11;  // # of intermeadiate profiles
  int opt_count = 0;
  double dAlpha = (1.0)/static_cast<double>(m-1);
  double t_fix = 0.5;
  const double advance_vec[] = {0., 0.5};
  std::vector<std::vector<double>> endData;
  std::vector<double> temp_end_data(6,55);

  // Wire parameters
  rbc::WireParams in_wire_params, out_wire_params, temp_wire_params;
  in_wire_params.H = 0.0;
  in_wire_params.V = 0.0;
  in_wire_params.M = 0.0;
  in_wire_params.BC = 0.0;
  in_wire_params.node = nNodes/2;
  
  // Wire cofnfiguration
  rbc::WireConfig in_wire_config, out_wire_config, temp_wire_config;
  in_wire_config.theta = std::vector<double>(nNodes, 0.);
  in_wire_config.xy = std::vector<double>(2*nNodes,0.);

  for(int j=0; j<nNodes; ++j)
    {
      in_wire_config.xy[2*j] = double(j)*static_cast<double>(blade_length/nElm);
      in_wire_config.xy[2*j] = redundentLen;
    }
  
  // looping over intermediate shapes
  double X[2]; int feasible = 0;
  std::vector<std::vector<double>> paramCont_interProfile;
  double box_min_y, elastica_max_y;
  std::vector<int> fix_node_details;
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
      fix_node_details.push_back(in_wire_params.node);
      
      // Advance the cut
      manager.AdvanceCut(optim_vars, in_wire_params, lbounds, ubounds, in_wire_config, intermediate_curve,
			 solve_params, nlopt_params, sdTOL, dtTOL, advance_vec, out_wire_params, out_wire_config, 0., 0.);
      
      std::cout << "\nis wire feasible : " <<  manager.IsWireFeasible(out_wire_config,sdTOL) << std::endl;
     
      // Output  = input for next cutting step
      in_wire_config.theta = out_wire_config.theta;
      in_wire_params = out_wire_params;

      PrintOpt(op_path+"opt-"+std::to_string(opt_count)+".dat", manager);
      CheckConstraint(target_curve,manager.GetBlade());
      
      //printig remaining area points
      PrintFile(op_path+"rem-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
      box_min_y = PrintBox(op_path+"box-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
      
      temp_end_data[0] = (out_wire_config.xy)[0];
      temp_end_data[1] = (out_wire_config.xy)[1];
      temp_end_data[2] = (out_wire_config.theta)[0];
      temp_end_data[3] = (out_wire_config.xy[2*(nNodes-1)]);
      temp_end_data[4] = (out_wire_config.xy[2*(nNodes-1)+1]);   
      temp_end_data[5] = (out_wire_config.theta[nNodes-1]);      
      endData.push_back(temp_end_data);
      ++opt_count;

      if(profile == m-3)
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
  rbc::WeightFunction wfunc(kval);
  int nFCLB = 1;
  
  int node_fix;
  double bounds[2] = {0, block_size};
  for(int nIter =0; nIter<nFCLB; ++nIter)
    {
      std::vector<double> uncut = uncutVec(manager.GetCutProfileMonitor());
      for(int k=0; k<(int)(uncut.size()/2); ++k)
	{
	  double u1 = uncut[2*k] - 0.05, u2 = uncut[2*k+1]+0.05;
	  if(u1 < 0 )
	    u1 = 0;
	  if(u2 > 1)
	    u2 = 1;
	  if(abs(u1 -u2) < 1e-3)
	    continue;
	  double weightRange[2] = {u1, u2};
	  t_fix = 0.5*(u1+ u2);
	  wfunc.SetIntervals(std::vector<std::pair<double,double>>{std::make_pair( u1, u2)});
	  weightFuncFixNode(nNodes,blade_length ,bounds, weightRange, init_curve,node_fix, t_fix);
	  temp_wire_params.node = node_fix;
	  target_curve.Evaluate(t_fix, temp_wire_params.xy);
	  fix_node_details.push_back(temp_wire_params.node);
	  
	  manager.AdvanceCut(optim_vars, temp_wire_params, lbounds, ubounds, temp_wire_config, target_curve,wfunc,			     
			     solve_params, nlopt_params, sdTOL, dtTOL, advance_vec, out_wire_params, out_wire_config, 0., 0.);     

	  PrintOpt(op_path+"opt-"+std::to_string(opt_count)+".dat", manager);
	  PrintFile(op_path+"rem-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
	  box_min_y = PrintBox(op_path+"box-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
	  ++opt_count;
	  PrintIntervals(manager.GetCutProfileMonitor());	  

	  temp_end_data[0] = (out_wire_config.xy)[0];
	  temp_end_data[1] = (out_wire_config.xy)[1];
	  temp_end_data[2] = (out_wire_config.theta)[0];
	  temp_end_data[3] = (out_wire_config.xy[2*nNodes-2]);
	  temp_end_data[4] = (out_wire_config.xy[2*nNodes-1]);   
	  temp_end_data[5] = (out_wire_config.theta[nNodes-1]);      
	  endData.push_back(temp_end_data);
	}
      uncut.clear();
      //temp_wire_config = out_wire_config;
      //temp_wire_params = out_wire_params;
    } 

  rbc::WireParams FCLB_in_wire_params, FCLB_out_wire_params;
  FCLB_in_wire_params.H = 0.0;
  FCLB_in_wire_params.V = 0.0;
  FCLB_in_wire_params.M = 0.0;
  FCLB_in_wire_params.BC = 0.0;
  
  // Wire cofnfiguration
  rbc::WireConfig FCLB_in_wire_config, FCLB_out_wire_config;
  FCLB_in_wire_config.theta = std::vector<double>(nNodes, 0.);
  FCLB_in_wire_config.xy = std::vector<double>(2*nNodes,0.);

  int nFCLB_WC = 2;
  
  for(int iter=0; iter<nFCLB_WC; ++iter)
    {
      std::vector<double> uncut = uncutVec(manager.GetCutProfileMonitor());
      for(int k=0; k<uncut.size()/2; ++k)
	{
	  double u1 = uncut[2*k] - 0.05, u2 = uncut[2*k+1]+0.05;
	  if(u1 < 0 )
	    u1 = 0;
	  if(u2 > 1)
	    u2 = 1;
	  if(abs(u1 -u2) < 1e-3)
	    continue;
	  double weightRange[2] = {u1, u2};
	  wfunc.SetIntervals(std::vector<std::pair<double,double>>{std::make_pair( u1, u2)});
	  weightFuncFixNode(nNodes,blade_length ,bounds, weightRange, init_curve,node_fix, t_fix);
	  FCLB_in_wire_params.node = node_fix;
	  t_fix = 0.5*(u1+ u2);

	  for(int j=0; j<nNodes; ++j)
	    {
	      FCLB_in_wire_config.xy[2*j] = double(j)*static_cast<double>(blade_length/nElm);
	      FCLB_in_wire_config.xy[2*j] = redundentLen;
	    }
      
	  for(int profile=0; profile<m; ++profile)
	    {
	      curve_combination.SetCombinationFactor(profile*dAlpha); // setting alpha
	      for(int j=0; j<nknots; ++j)	
		{
		  curve_combination.Evaluate(tvec[j],X);
		  inter_xy[2*j] = X[0];
		  inter_xy[2*j+1] = X[1];
		}
	      intermediate_curve.SetControlPoints(inter_xy);
	      intermediate_curve.Evaluate(t_fix, FCLB_in_wire_params.xy);
	  
	      // Advance the cut
	      manager.AdvanceCut(optim_vars, FCLB_in_wire_params, lbounds, ubounds, FCLB_in_wire_config, intermediate_curve,wfunc,
				 solve_params, nlopt_params, sdTOL, dtTOL, advance_vec, FCLB_out_wire_params, FCLB_out_wire_config, 0., 0.);      
     
	      // Output  = input for next cutting step
	      FCLB_in_wire_config.theta = FCLB_out_wire_config.theta;
	      FCLB_in_wire_params = FCLB_out_wire_params;

	      CheckConstraint(target_curve,manager.GetBlade());
      
	      //printig remaining area points
	      if(profile == m-1)
		{
		  PrintOpt(op_path+"opt-"+std::to_string(opt_count)+".dat", manager);
		  PrintFile(op_path+"rem-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
		  box_min_y = PrintBox(op_path+"box-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
		  ++opt_count;
		  fix_node_details.push_back(FCLB_in_wire_params.node);
		  PrintIntervals(manager.GetCutProfileMonitor());	  		  
		}
	    }
	}
    }
  std::vector<std::vector<double>> op1 = xy_map_alphaR(endData, blade_length, redundentLen);

  std::fstream data;
  data.open("../Data/fix_node.dat",std::ios::out);
  for(int i=0; i<fix_node_details.size(); ++i)
    data << fix_node_details[i] << "\n";
  data.close();
}// main function ends

//**************************Functions***************************//
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


void arrangeProfile(const std::vector<std::vector<double>> op1, std::vector<std::vector<double>>& arrangedData)
{
  // 4 possibility
  bool end_0, end_1; // bool =0 >>> -ve, 
}
