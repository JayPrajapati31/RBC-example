// Sriramajayam

#ifndef JSON_READER_H
#define JSON_READER_H

#include "json.hpp"
#include <iostream>

namespace JSON_Read
{
  struct json_BladeParams
  {
    const std::string tag = "blade params";
    double length;
    double EI;
    int nodes;
    friend std::istream& operator >> (std::istream &in, json_BladeParams& params);
  };

  struct json_NLoptParams
  {
    const std::string tag = "nlopt params";
    double relative_ftol;
    double relative_xtol;
    friend std::istream& operator >> (std::istream &in, json_NLoptParams& params);
  };

  struct json_IPoptParams
  {
    const std::string tag = "ipopt params";
    std::string linear_solver;
    double obj_scale_factor;
    std::string hessian_approx;
    std::string warm_start_init;
    double warm_start_bound;  
    double warm_start_mult_bound;  
    friend std::istream& operator >> (std::istream &in, json_IPoptParams& params);
  };

  struct json_TargetShape
  {
    const std::string tag = "target shape";
    std::string filename;
    int nSamples_per_interval;
    friend std::istream& operator >> (std::istream &in, json_TargetShape& params);
  };

  struct json_Block
  {
    const std::string tag = "virtual block";
    double center[2];
    double size;
    friend std::istream& operator >> (std::istream &in, json_Block& params);
  };

  struct json_CutTol
  {
    const std::string tag = "cut tolerances";
    double sdTol;
    int nElm_blade;
    friend std::istream& operator >> (std::istream &in, json_CutTol& params);
  };

  struct json_ElasticaSolverParams
  {
    const std::string tag = "elastica solver params";
    double tol;
    double rScale;
    double dofScale;
    int maxIter;
    bool verbosity;
    friend std::istream& operator >> (std::istream &in, json_ElasticaSolverParams& params);
  };

  struct json_OptParams
  {
    const std::string tag = "optimization parameters";
    bool H;
    bool V;
    bool M;
    bool BC;
    bool axial;
    bool radial;
    friend std::istream& operator >> (std::istream &in, json_OptParams& params);
  };

  struct json_OptParamBounds
  {
    const std::string tag = "parameter bounds";
    double Hbound[2];
    double Vbound[2];
    double Mbound[2];
    double BCbound[2];
    double axialBound[2];
    double radialBound[2];
    friend std::istream& operator >> (std::istream &in, json_OptParamBounds& params);
  };

  struct json_RoughCutParams
  {
    const std::string tag = "rough cut";
    int nProfile;
    double t_fix; // on profile
    double s_fix; // on blade
    friend std::istream& operator >> (std::istream &in, json_RoughCutParams& params);
  };

  struct json_InitGuess
  {
    const std::string tag = "initial params";
    double H;
    double V;
    double M;
    double BC;
    friend std::istream& operator >> (std::istream &in, json_InitGuess& params);
  };

  struct json_WfuncParams
  {
    const std::string tag = "weight function param";
    double kval;
    friend std::istream& operator >> (std::istream &in, json_WfuncParams& params);
  };

}
//void ReadJSON(const std::string filename);

#endif
