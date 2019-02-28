#ifndef _HARMONIC_SOLVER_
#define _HARMONIC_SOLVER_

#include "head.h"
#include "myvector.h"
#include "object.h"

class harmonic_solver
{
 public:
  double dmetric;
  double scale;//force need to be scale
  harmonic_solver(){}
  ~harmonic_solver(){}
  harmonic_solver(boost::property_tree::ptree &para_tree);
  int default_solve(boost::property_tree::ptree &para_tree);
  int solve(boost::property_tree::ptree &para_tree);
  int configMaterial(object *myobject,const double (&material_para)[2]); 
  int configForce(object *myobject, myvector e1, myvector e2);
};
#endif
