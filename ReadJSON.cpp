// Sriramajayam


#include "ReadJSON.h"
#include <fstream>
#include <cassert>
#include <iostream>

namespace JSON_Read
{
  std::istream& operator >> (std::istream &in, json_BladeParams& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("length"));
    assert(J[params.tag].contains("EI"));
    assert(J[params.tag].contains("nodes"));

    // read
    J[params.tag]["length"].get_to(params.length);
    J[params.tag]["EI"].get_to(params.EI);
    J[params.tag]["nodes"].get_to(params.nodes);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_NLoptParams& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("relative ftol"));
    assert(J[params.tag].contains("relative xtol"));

    // read
    J[params.tag]["relative ftol"].get_to(params.relative_ftol);
    J[params.tag]["relative xtol"].get_to(params.relative_xtol);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_IPoptParams& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("linear_solver"));
    assert(J[params.tag].contains("objective scale factor"));
    assert(J[params.tag].contains("hessian approximation"));
    assert(J[params.tag].contains("warm start initial point"));
    assert(J[params.tag].contains("warm start bound"));
    assert(J[params.tag].contains("warm start multiplier bound"));

    // read
    J[params.tag]["linear_solver"].get_to(params.linear_solver);
    J[params.tag]["objective scale factor"].get_to(params.obj_scale_factor);
    J[params.tag]["hessian approximation"].get_to(params.hessian_approx);
    J[params.tag]["warm start initial point"].get_to(params.warm_start_init);
    J[params.tag]["warm start bound"].get_to(params.warm_start_bound);
    J[params.tag]["warm start multiplier bound"].get_to(params.warm_start_mult_bound);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_TargetShape& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("filename"));
    assert(J[params.tag].contains("num samples per interval"));

    // read
    J[params.tag]["filename"].get_to(params.filename);
    J[params.tag]["num samples per interval"].get_to(params.nSamples_per_interval);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_Block& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("center"));
    assert(J[params.tag].contains("size"));

    // read
    J[params.tag]["center"].get_to(params.center);
    J[params.tag]["size"].get_to(params.size);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_CutTol& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("distance tolerance"));
    assert(J[params.tag].contains("elements on blade"));

    // read
    J[params.tag]["distance tolerance"].get_to(params.sdTol);
    J[params.tag]["elements on blade"].get_to(params.nElm_blade);

    return in;
  }


  std::istream& operator >> (std::istream &in, json_OptParams& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("H"));
    assert(J[params.tag].contains("V"));
    assert(J[params.tag].contains("M"));
    assert(J[params.tag].contains("BC"));
    assert(J[params.tag].contains("axial"));
    assert(J[params.tag].contains("radial"));
  
    // read
    J[params.tag]["H"].get_to(params.H);
    J[params.tag]["V"].get_to(params.V);
    J[params.tag]["M"].get_to(params.M);
    J[params.tag]["BC"].get_to(params.BC);
    J[params.tag]["axial"].get_to(params.axial);
    J[params.tag]["radial"].get_to(params.radial);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_OptParamBounds& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("H"));
    assert(J[params.tag].contains("V"));
    assert(J[params.tag].contains("M"));
    assert(J[params.tag].contains("BC"));
    assert(J[params.tag].contains("axial"));
    assert(J[params.tag].contains("radial"));
  
    // read
    J[params.tag]["H"].get_to(params.Hbound);
    J[params.tag]["V"].get_to(params.Vbound);
    J[params.tag]["M"].get_to(params.Mbound);
    J[params.tag]["BC"].get_to(params.BCbound);
    J[params.tag]["axial"].get_to(params.axialBound);
    J[params.tag]["radial"].get_to(params.radialBound);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_RoughCutParams& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("number of profiles"));
    assert(J[params.tag].contains("mapped parameter on profile"));
    assert(J[params.tag].contains("mapped node on blade"));

    // read
    J[params.tag]["number of profiles"].get_to(params.nProfile);
    J[params.tag]["mapped parameter on profile"].get_to(params.t_fix);
    J[params.tag]["mapped node on blade"].get_to(params.s_fix);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_InitGuess& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("H"));
    assert(J[params.tag].contains("V"));
    assert(J[params.tag].contains("M"));
    assert(J[params.tag].contains("BC"));

    // read
    J[params.tag]["H"].get_to(params.H);
    J[params.tag]["V"].get_to(params.V);
    J[params.tag]["M"].get_to(params.M);
    J[params.tag]["BC"].get_to(params.BC);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_WfuncParams& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    //assert(J[params.tag].contains("weight function param"));


    // read
    J[params.tag].get_to(params.kval);

    return in;
  }

  std::istream& operator >> (std::istream &in, json_ElasticaSolverParams& params)
  {
    assert(params.tag.empty()==false);

    // rewind to the top
    in.seekg(0, std::ios::beg);

    // read until the component is found
    nlohmann::json J;
    while(J.empty() && in.good())
      {
	in >> J;
	if(J.contains(params.tag)==false)
	  J.clear();
      }

    // check fields
    assert(J.empty()==false);
    assert(J[params.tag].contains("convergence tol"));
    assert(J[params.tag].contains("residual scaling"));
    assert(J[params.tag].contains("dof scaling"));
    assert(J[params.tag].contains("maximum iterations"));
    assert(J[params.tag].contains("verbosity"));

    // read
    //J[params.tag]["convergence tol"].get_to(params.tol);
    J[params.tag]["residual scaling"].get_to(params.rScale);
    J[params.tag]["dof scaling"].get_to(params.dofScale);
    J[params.tag]["maximum iterations"].get_to(params.maxIter);
    J[params.tag]["verbosity"].get_to(params.verbosity);

    return in;
  }

  /*void ReadJSON(const std::string filename)
    {
    std::fstream stream;
    stream.open(filename.c_str(), std::ios::in);
    assert(stream.good());

    json_BladeParams blade_params;
    json_NLoptParams nlopt_params;
    json_IPoptParams ipopt_params;
    json_TargetShape targetShape;
    json_Block block;
    json_CutTol cutTol;
    json_ElasticaSolverParams elasticaSolverParams;
    json_OptParams optParams;
    json_OptParamBounds optParamBounds;
    json_RoughCutParams roughCutParams;
    json_InitGuess initGuess;
    json_WfuncParams wfuncParams;
    stream >> blade_params;
    stream >> nlopt_params;
    stream >> ipopt_params;
    stream >> targetShape;
    stream >> block;
    stream >> cutTol;
    stream >> optParams;
    stream >> optParamBounds;
    stream >> roughCutParams;
    stream >> initGuess;
    stream >> wfuncParams;
    stream >> elasticaSolverParams;

    std::cout<<"\nBlade parameters: "<<blade_params.length<<", "<<blade_params.EI<<", "<<blade_params.nodes<<std::endl;

    std::cout<<"\nNLopt Parameters: "<<nlopt_params.relative_ftol<<", "<<nlopt_params.relative_xtol<<std::endl;

    std::cout<<"\nIPopt Parameters: "<< ipopt_params.linear_solver<<", "<<ipopt_params.obj_scale_factor<<", "<< ipopt_params.hessian_approx << ", "<< ipopt_params.warm_start_init<< ", "<< ipopt_params.warm_start_bound << ", "<< ipopt_params.warm_start_mult_bound << std::endl;

    std::cout<<"\nTargte shape params: "<<targetShape.filename<<", "<<targetShape.nSamples_per_interval<<std::endl;
  
    std::cout << "\nBloxk params" <<block.center[0] << ", " << block.center[1] << ", " << block.size << std::endl;

    std::cout<<"\ncutting tolerances: "<<cutTol.sdTol<<", "<<cutTol.nElm_blade<<std::endl;

    std::cout << "\nOpt params selected " << optParams.H << ", " << optParams.V << ", " << optParams.M << ", " << optParams.BC << ", " << optParams.axial << ", " << optParams.radial << std::endl;
    stream.close();

    std::cout << "\nOpt param bounds " << "[" << optParamBounds.Hbound[0] << ", " << optParamBounds.Hbound[1] << "]\n"
    << "[" << optParamBounds.Vbound[0] << ", " << optParamBounds.Vbound[1] << "]\n"
    << "[" << optParamBounds.Mbound[0] << ", " << optParamBounds.Mbound[1] << "]\n"
    << "[" << optParamBounds.BCbound[0] << ", " << optParamBounds.BCbound[1] << "]\n"
    << "[" << optParamBounds.axialBound[0] << ", " << optParamBounds.axialBound[1] << "]\n"
    << "[" << optParamBounds.radialBound[0] << ", " << optParamBounds.radialBound[1] << "]"
    << std::endl;

    std::cout << "\nRough cut params " << roughCutParams.nProfile << ", " << roughCutParams.t_fix << ", " << roughCutParams.s_fix << std::endl;

    std::cout << "\nInit value of opt params " << initGuess.H << ", \n"
    << initGuess.V << ", \n"
    << initGuess.M << ", \n"
    << initGuess.BC << std::endl;

    std::cout << "\nWfunc params " << wfuncParams.kval << std::endl;  

    std::cout << "\elastica solver params" << elasticaSolverParams.tol << ", " << elasticaSolverParams.rScale << ", " << elasticaSolverParams.dofScale << ", " << elasticaSolverParams.maxIter << ", " << elasticaSolverParams.verbosity << std::endl;
  
    }*/

}
