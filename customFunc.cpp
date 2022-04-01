#include "customFunc.h"
#include <iostream>
#include <fstream>

bool CheckConstraint(rbc::ParametricCurve& splineObj ,const rbc::Elastica& elasticaObj)
{
  int n = elasticaObj.GetNumNodes();
  std::vector<double> xy_sol(2*n);
  xy_sol =  elasticaObj.GetCartesianCoordinates();
  std::vector<double> op(n);
  double coords[2], sd;
  bool flag;
  
  for(int i=0; i<n; ++i)
    {
      coords[0] = xy_sol[2*i];
      coords[1] = xy_sol[2*i+1];
      splineObj.GetSignedDistance(coords, sd, nullptr);
      op[i] = sd;
      if(sd>=0)
	{
	  flag = false;
	  std::cout << " \n Node " << i  << " sd : " << sd;
	  break;
	}
    }
  return flag;
} 


void PrintIntervals(const rbc::CutProfileMonitor& monitor)
{
  auto& cut_intervals = monitor.GetCutIntervals();
  auto& uncut_intervals = monitor.GetUncutIntervals();

  std::cout<<"\nUncut intervals: ";
  for(auto& it:uncut_intervals)
    std::cout<<"("<<it.lower()<<","<<it.upper()<<") "<<std::endl;
}


std::vector<double> uncutVec(const rbc::CutProfileMonitor& monitor)
{
  std::vector<double> uncut;
  auto& uncut_intervals = monitor.GetUncutIntervals();
  for(auto& it:uncut_intervals)
    {
      uncut.push_back( it.lower() );
      uncut.push_back( it.upper() );
    }
  
  return uncut; 
}

/*template<class T>
void PrintFile(std::string filename, T& obj)
{
  std::fstream data;
  data.open(filename, std::ios::out);
  data << obj;
  data.close();
}*/

void PrintOpt(std::string filename, const rbc::ProfileCuttingManager& manager)
{
  std::fstream data;
  data.open(filename, std::ios::out);
  std::vector<double> xy = manager.GetBlade().GetCartesianCoordinates();
  std::vector<double> coords = manager.GetBlade().GetCoordinates();
  std::vector<double> states = manager.GetBlade().GetStateField();
  for(int k=0; k<manager.GetBlade().GetNumNodes(); ++k)
    {
      data << coords[k]<<"\t" << states[k] << "\t" << xy[2*k] << "\t" << xy[2*k+1] << std::endl;
    }
  data.close();
}


double PrintBox(std::string fname, const rbc::CutAreaMonitor& monitor)
{
  boost_box bounding_box;
  boost::geometry::envelope(monitor.GetUncutPolygon(), bounding_box);
  double min_x = bounding_box.min_corner().get<0>();
  double min_y = bounding_box.min_corner().get<1>();
  double max_x = bounding_box.max_corner().get<0>();
  double max_y = bounding_box.max_corner().get<1>();
  std::vector<double> box = {min_x, min_y,
			     min_x, max_y,
			     max_x, max_y,
			     max_x, min_y,
			     min_x, min_y};
  assert(box.size() == 10);
  std::fstream data;
  data.open(fname, std::ios::out);
  for(int i=0; i<box.size()/2; ++i)
    data << box[2*i] << "\t" << box[2*i+1] << std::endl;
  data.close();
  return min_y;
}

// takes xy data of the ends and returns rad, alpha and theta
std::vector<std::vector<double>> xy_map_alphaR(std::vector<std::vector<double>> end_data, double bladeLen, double redundent_length)
{
  std::vector<std::vector<double>> op(end_data.size(), std::vector<double>(6,0));
  std::vector<std::vector<double>> op1(end_data.size(), std::vector<double>(6,0));

  double x0_t, y0_t, x0_b, y0_b;
  double dx_t, dy_t, dx_b, dy_b , ratio_t, ratio_b;
  x0_t = end_data[0][0];
  y0_t = end_data[0][1]-redundent_length;
  x0_b = end_data[0][3];
  y0_b = end_data[0][4]-redundent_length;

  //std::cout << x0_t << "\t" << y0_t <<"\t" << x0_b << "\t" << y0_b << std::endl;
  
  for(int i=0; i< end_data.size(); ++i)
    {
      dx_t = end_data[i][0] - x0_t;
      dy_t = end_data[i][1] - y0_t;
      ratio_t = (dy_t/ dx_t);
      if(dx_t == 0)
	ratio_t = HUGE_VAL;
      
      //radius
      op[i][0] = pow( pow(dx_t,2) + pow(dy_t,2) ,0.5); // radius
      op1[i][0] = pow( pow(dx_t,2) + pow(dy_t,2) ,0.5); // radius

      // alpha
      op[i][1] = 0.5*M_PI - atan2(dy_t, dx_t); // for experiments
      op1[i][1] = atan2(dy_t, dx_t); // for plotting

      //theta
      op[i][2] = end_data[i][2]; // theta
      op1[i][2] = op[i][2];

      
      dx_b = end_data[i][3] - x0_b;
      dy_b = end_data[i][4] - y0_b;
      ratio_b = dy_b/dx_b;
      if(dx_b == 0)
	ratio_b = HUGE_VAL;

      //radius
      op[i][3] = pow(pow(dx_b,2) + pow(dy_b,2) ,0.5); // radius
      op1[i][3] = pow(pow(dx_b,2) + pow(dy_b,2) ,0.5); // radius

      //alpha
      op[i][4] =  atan2(dy_b, dx_b) - 0.5*M_PI; // alpha
      op1[i][4] = atan2(dy_b, dx_b);

      //theta
      op[i][5] = end_data[i][5]; // theta
      op1[i][5] = op[i][5];
    }
   
  std::fstream data;
  data.open("../Data/r_alpha_theta.dat", std::ios::out);
  for(int i=0; i<end_data.size(); ++i)
    {
      data << "{";
      for(int j =0; j< 6; ++j)
	{
	  data << op[i][j] << "\t";
	  if(j<5)
	    data << ",";
	}
      data << "},"<<   std::endl;
    }
  data.close();

  data.open("../Data/exp.dat", std::ios::out);
  for(int i=0; i<end_data.size(); ++i)
    {
      for(int j =0; j< 6; ++j)
	{
	  data << op[i][j] << "\t";
	}
      data << std::endl;
    }
  data.close();
  
  data.open("../Data/sim.dat", std::ios::out);
  for(int i=0; i<end_data.size(); ++i)
    {
      for(int j =0; j< 6; ++j)
	{
	  data << op1[i][j] << "\t";
	}
      data << std::endl;
    }
  data.close();
  
  return op1;
}

int ElasticaMax(const rbc::Elastica& elastica, double* X)
{
  std::vector<double> xy = elastica.GetCartesianCoordinates();
  std::vector<double> y;
  
  for(int i=0; i<elastica.GetNumNodes(); ++i)
    y.push_back(xy[2*i+1]);

  int node = max_element(y.begin(), y.end()) - y.begin();
  X[0] = xy[2*node];
  X[1] = xy[2*node+1];
  // std::cout << "\n" << X[0] << "\t"<< X[1] << std::endl;
  return node;
}

// intermediate profile continuation
int intermediateProfContinuation(const std::vector<std::vector<double>> paramVec, int upto_profile,const double blade_length, int counter)
{
  const int nNodes = 101;
  double h = blade_length/static_cast<double>(nNodes-1);
  std::vector<double> coords(nNodes);
  for(int i=0; i<nNodes; ++i)
    {    
      coords[i] = i*h;
    }
  rbc::LoadParams load_params({.HVal=0., .VVal=0., .MDof=nNodes-1, .MVal=2.*M_PI});  
  std::pair<int,double> bc_params(0, 0.);
  rbc::OriginParams origin_params{.node=0, .coords={0.,0.}};
  rbc::SolveParams solve_params({.EPS = 1.e-10, .resScale = 1.,	.dofScale = 1.,	.nMaxIters = 25, .verbose = false});
  const double EI = 1.;
  rbc::Elastica elastica(coords, EI);
  
  int n = paramVec.size();
  assert(n!=0);

  int ContSteps = 5;
  double dLambda = (1.0)/static_cast<double>(ContSteps-1);
  double lambda = 0.;

  // do continuation till "upto_profile" intermediate profile
  for(int i=0; i< upto_profile; ++i)
      for(int j=0; j<ContSteps; ++j)
	{
	  lambda = j*dLambda;
	  load_params.HVal = lambda*paramVec[i+1][0] + (1.0-lambda)*paramVec[i][0];
	  load_params.VVal = lambda*paramVec[i+1][1] + (1.0-lambda)*paramVec[i][1];
	  load_params.MVal = lambda*paramVec[i+1][2] + (1.0-lambda)*paramVec[i][2];
	  bc_params.second = lambda*paramVec[i+1][3] + (1.0-lambda)*paramVec[i][3];
	  origin_params.node =(int) paramVec[0][4];
	  origin_params.coords[0] = lambda*paramVec[i+1][5] + (1.0-lambda)*paramVec[i][5];
	  origin_params.coords[1] = lambda*paramVec[i+1][6] + (1.0-lambda)*paramVec[i][6];
	  elastica.ComputeState(load_params,bc_params,origin_params,solve_params);

	  PrintFile("../Continuation Data/cont-"+std::to_string(counter)+".dat", elastica);
	  ++counter;
	}
  
  return counter;
}

void readTarget(std::string filename, std::vector<double>& xy, std::vector<double>& tvec)
{
  double temp[3];
  std::fstream data;
  data.open(filename,std::ios::in);
  std::string ln;
  while(std::getline(data,ln,'\n'))
    {
      std::istringstream iss(ln);
      std::string field;
      int j = 0;
      while(getline(iss,field,'\t'))
	{
	  temp[j] = std::stod(field); ++j;
	}
      tvec.push_back(temp[0]);
      xy.push_back(temp[1]);
      xy.push_back(temp[2]);
    }
  data.close();
}
