#ifndef CUSTOM_FUNC_H
#define CUSTOM_FUNC_H

#include <rbc_ParametricCurve.h>
#include <rbc_CubicBezierCurve.h>
#include <rbc_CubicSplineCurve.h>
#include <rbc_Elastica.h>
#include <rbc_ProfileCuttingManager.h>
#include <rbc_ProfileCuttingManager_impl.h>
#include <rbc_CutAreaMonitor.h>
#include <rbc_CutProfileMonitor.h>
#include <vector>
#include <fstream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/geometries.hpp>
typedef boost::geometry::model::box<rbc::boost_point2D> boost_box;

template<class T>
void PrintFile(std::string filename, T& obj)
{
  std::fstream data;
  data.open(filename, std::ios::out);
  data << obj;
  data.close();
}

bool CheckConstraint(rbc::ParametricCurve& splineObj ,const rbc::Elastica& elasticaObj);

void PrintIntervals(const rbc::CutProfileMonitor& monitor);

std::vector<double> uncutVec(const rbc::CutProfileMonitor& monitor);

void PrintOpt(std::string filename, const rbc::ProfileCuttingManager& manager);

double PrintBox(std::string fname, const rbc::CutAreaMonitor& monitor);
 
std::vector<std::vector<double>> xy_map_alphaR(std::vector<std::vector<double>> end_data, double bladeLen, double redundent_length);

int ElasticaMax(const rbc::Elastica& elastica, double* X);

int intermediateProfContinuation(const std::vector<std::vector<double>> paramVec, int upto_profile,const double blade_length, int counter);

void readTarget(std::string fileName, std::vector<double>& xy, std::vector<double>& TVEC);

#endif
