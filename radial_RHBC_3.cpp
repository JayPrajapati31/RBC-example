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
#include "ReadJSON.h"

JSON_Read::json_BladeParams json_blade_params;
JSON_Read::json_NLoptParams json_nlopt_params;
JSON_Read::json_IPoptParams json_ipopt_params;
JSON_Read::json_TargetShape json_targetShape;
JSON_Read::json_Block json_block;
JSON_Read::json_CutTol json_cutTol;
JSON_Read::json_ElasticaSolverParams json_elasticaSolverParams;
JSON_Read::json_OptParams json_optParams;
JSON_Read::json_OptParamBounds json_optParamBounds;
JSON_Read::json_RoughCutParams json_roughCutParams;
JSON_Read::json_InitGuess json_initGuess;
JSON_Read::json_WfuncParams json_wfuncParams;

void ReadJSON(const std::string filename);
void weightFuncFixNode(const int nNodes,const double blade_len,const double* bounds, const double* weightRange, const rbc::ParametricCurve& init,  int& node, double& tparam);

int main(int argc, char** argv)
{
  ReadJSON("../input.json");

  // Initialize PETSc
  PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

  // Create a uniform 1D mesh on the interval [0,1]
  const double kval = json_wfuncParams.kval;
  double blade_length = json_blade_params.length;
  const int nNodes = json_blade_params.nodes;
  const double redundentLen = 0.2525;
  const int nElm = nNodes-1;

  // Target shape
  rbc::NLopt_Params nlopt_params;
  nlopt_params.set_ftol_abs = false ;nlopt_params.set_ftol_rel = true;
  nlopt_params.set_xtol_abs = false; nlopt_params.set_xtol_rel = true;
  const double ftol = json_nlopt_params.relative_ftol;
  const double xtol = json_nlopt_params.relative_xtol;
  nlopt_params.ftol_abs = ftol;
  nlopt_params.ftol_rel = ftol;
  nlopt_params.xtol_abs = xtol;
  nlopt_params.xtol_rel = xtol;

  rbc::Ipopt_Params ipopt_params; 

  std::vector<double> tvec, xy;
  readTarget(json_targetShape.filename, xy, tvec);
  assert(tvec.size() == xy.size()/2);
  int nknots = tvec.size();
  int nSample = json_targetShape.nSamples_per_interval;
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
  
  // linear combination of profiles
  rbc::ParametricCurveCombination curve_combination(target_curve,init_curve,100,nlopt_params);

  // Profile cutting manager
  rbc::ProfileCuttingManager manager(target_curve, std::vector<double>{0.,redundentLen,
									 0., y_uBound + 0.1,
									 block_size, y_uBound+0.1,
									 block_size, redundentLen});
  manager.CreateBlade(blade_length, nElm);  
    
  // Constants
  const double sdTOL = json_cutTol.sdTol;  // determines which points on the target has been cut
  const double dtTOL = 0.005; // size of smallest cut = size of 4 element
    
  // Solve parameters
  rbc::SolveParams solve_params({.EPS = json_elasticaSolverParams.tol, .resScale = json_elasticaSolverParams.rScale, .dofScale = json_elasticaSolverParams.dofScale, .nMaxIters = json_elasticaSolverParams.maxIter, .verbose = json_elasticaSolverParams.verbosity});
  
  // define optimization variables
  std::array<bool,6> optim_vars{json_optParams.H, json_optParams.V, json_optParams.M, json_optParams.BC, json_optParams.axial, json_optParams.radial}; // H, V, M, BC, r=y
  
  // bounds for optimization parameters
  std::vector<double> lbounds(6), ubounds(6),ip_lbounds(5), ip_ubounds(5) ;
  lbounds[0] = json_optParamBounds.Hbound[0];
  lbounds[1] = json_optParamBounds.Vbound[0];
  lbounds[2] = json_optParamBounds.Mbound[0];
  lbounds[3] = json_optParamBounds.BCbound[0];
  lbounds[4] = json_optParamBounds.radialBound[0];

  ubounds[0] = json_optParamBounds.Hbound[1];
  ubounds[1] = json_optParamBounds.Vbound[1];
  ubounds[2] = json_optParamBounds.Mbound[1];
  ubounds[3] = json_optParamBounds.BCbound[1];
  ubounds[4] = json_optParamBounds.radialBound[1];
  
  // path to store results
  std::string op_path = "/home/jay/CodeRepos/rbc/examples/axisymm/Radial RHBC/Data/"; 
  
  // alpha and deltaAlpha
  int m = json_roughCutParams.nProfile;  // # of intermeadiate profiles
  int opt_count = 0;
  double dAlpha = (1.0)/static_cast<double>(m-1);
  double t_fix = json_roughCutParams.t_fix;
  const double advance_vec[] = {0., 0.5};
  std::vector<std::vector<double>> endData;
  std::vector<double> temp_end_data(6,55);

  // Wire parameters
  rbc::WireParams in_wire_params, out_wire_params, temp_wire_params;
  in_wire_params.H = json_initGuess.H;
  in_wire_params.V = json_initGuess.V;
  in_wire_params.M = json_initGuess.M;
  in_wire_params.BC = json_initGuess.BC;
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
  std::vector<double> fix_node_details;
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
      
      // Advance the cut
      manager.AdvanceCut(optim_vars, in_wire_params, lbounds, ubounds, in_wire_config, intermediate_curve,
			 solve_params, nlopt_params, sdTOL, dtTOL, advance_vec, out_wire_params, out_wire_config, 0., 0.);
      
      std::cout << "\nis wire feasible : " <<  manager.IsWireFeasible(out_wire_config,sdTOL) << std::endl;
     
      // Output  = input for next cutting step
      in_wire_config.theta = out_wire_config.theta;
      in_wire_params = out_wire_params;
      fix_node_details.push_back(in_wire_params.node);
      
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

      if(profile == m-2)
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
  int nFCLB = 2;
  
  /*int node_fix;
  double bounds[2] = {0, block_size};
  for(int nIter =0; nIter<nFCLB; ++nIter)
    {
      std::vector<double> uncut = uncutVec(manager.GetCutProfileMonitor());
      for(int k=0; k<(int)(uncut.size()/2); ++k)
	{       
	  double weightRange[2] = {uncut[2*k], uncut[2*k+1]};
	  t_fix = 0.5*(uncut[2*k]+ uncut[2*k+1]);
	  wfunc.SetIntervals(std::vector<std::pair<double,double>>{std::make_pair( uncut[2*k], uncut[2*k+1])});
	  weightFuncFixNode(nNodes,blade_length ,bounds, weightRange, init_curve,node_fix, t_fix);
	  temp_wire_params.node = node_fix;
	  target_curve.Evaluate(t_fix, temp_wire_params.xy);

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
	  fix_node_details.push_back(out_wire_params.node);
	}
      uncut.clear();
      temp_wire_config = out_wire_config;
      temp_wire_params = out_wire_params;
    }*/


  /***************************FCLB*******************************/
  rbc::WireParams FCLB_in_wire_params, FCLB_out_wire_params;
    FCLB_in_wire_params.H = json_initGuess.H;
    FCLB_in_wire_params.V = json_initGuess.V;
    FCLB_in_wire_params.M = json_initGuess.M;
    FCLB_in_wire_params.BC = json_initGuess.BC;
    FCLB_in_wire_params.node = nNodes/2;
  
    // Wire cofnfiguration
    rbc::WireConfig FCLB_in_wire_config, FCLB_out_wire_config;
    FCLB_in_wire_config.theta = std::vector<double>(nNodes, 0.);
    FCLB_in_wire_config.xy = std::vector<double>(2*nNodes,0.);
    

     double bounds[2] = {0, block_size};
  for(int nIter =0; nIter<nFCLB; ++nIter)
    {
      for(int j=0; j<nNodes; ++j)
      {
	FCLB_in_wire_config.xy[2*j] = double(j)*static_cast<double>(blade_length/nElm);
	FCLB_in_wire_config.xy[2*j] = redundentLen;
      }
      std::vector<double> uncut = uncutVec(manager.GetCutProfileMonitor());
      for(int k=0; k<(int)(uncut.size()/2); ++k)
	{
	  if(uncut[2*k]==uncut[2*k+1])
	    continue;
	  for(int profile=0;profile<m; ++profile)
	    {
	      curve_combination.SetCombinationFactor(profile*dAlpha); // setting alpha
	      for(int j=0; j<nknots; ++j)	
		{
		  curve_combination.Evaluate(tvec[j],X);
		  inter_xy[2*j] = X[0];
		  inter_xy[2*j+1] = X[1];
		}
	      intermediate_curve.SetControlPoints(inter_xy);	      
	      
	      double weightRange[2] = {uncut[2*k], uncut[2*k+1]};
	      wfunc.SetIntervals(std::vector<std::pair<double,double>>{std::make_pair( uncut[2*k], uncut[2*k+1])});
	      weightFuncFixNode(nNodes,blade_length ,bounds, weightRange, init_curve,FCLB_in_wire_params.node, t_fix);
	      t_fix = 0.5*(uncut[2*k]+ uncut[2*k+1]);
	      target_curve.Evaluate(t_fix, FCLB_in_wire_params.xy);

	      manager.AdvanceCut(optim_vars, FCLB_in_wire_params, lbounds, ubounds, FCLB_in_wire_config, target_curve,wfunc,			     
				 solve_params, nlopt_params, sdTOL, dtTOL, advance_vec, FCLB_out_wire_params, FCLB_out_wire_config, 0., 0.);	  
	  
	      PrintOpt(op_path+"opt-"+std::to_string(opt_count)+".dat", manager);
	      PrintFile(op_path+"rem-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
	      box_min_y = PrintBox(op_path+"box-"+std::to_string(opt_count)+".dat", manager.GetCutAreaMonitor());
	      ++opt_count;
	      PrintIntervals(manager.GetCutProfileMonitor());	
	  
	      temp_end_data[0] = (FCLB_out_wire_config.xy)[0];
	      temp_end_data[1] = (FCLB_out_wire_config.xy)[1];
	      temp_end_data[2] = (FCLB_out_wire_config.theta)[0];
	      temp_end_data[3] = (FCLB_out_wire_config.xy[2*nNodes-2]);
	      temp_end_data[4] = (FCLB_out_wire_config.xy[2*nNodes-1]);   
	      temp_end_data[5] = (FCLB_out_wire_config.theta[nNodes-1]);	  
	      endData.push_back(temp_end_data);
	     
	      FCLB_in_wire_params = FCLB_out_wire_params;
	      fix_node_details.push_back(FCLB_in_wire_params.node);
	    }
	}
      uncut.clear();
    }

  std::fstream data;
  data.open("../Data/fix_node.dat", std::ios::out);
  for(int i=0; i<fix_node_details.size(); ++i)
    data << fix_node_details[i] << "\n";
  data.close();
  std::vector<std::vector<double>> op1 = xy_map_alphaR(endData, blade_length, redundentLen);
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
  double x_mid = t_mid*(bounds[1]-bounds[0]) + bounds[0];
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

void ReadJSON(const std::string filename)
{
  std::fstream stream;
  stream.open(filename.c_str(), std::ios::in);
  assert(stream.good());

  stream >> json_blade_params;
  stream >> json_nlopt_params;
  stream >> json_ipopt_params;
  stream >> json_targetShape;
  stream >> json_block;
  stream >> json_cutTol;
  stream >> json_optParams;
  stream >> json_optParamBounds;
  stream >> json_roughCutParams;
  stream >> json_initGuess;
  stream >> json_wfuncParams;
  stream >> json_elasticaSolverParams;
}
