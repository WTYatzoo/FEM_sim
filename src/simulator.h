#ifndef _SIMULATOR_
#define _SIMULATOR_

#include "head.h"
#include "object.h"
class simulator
{
 public:
  std::string kind; //default for test or use data from file 
  std::string simulation_type;
  size_t frame;
  double dt;
  std::string force_function;
  double gravity;
  double density;

  int line_search;
  double weight_line_search;
  std::string constitutive_model;
  std::string out_dir;
  simulator(){}
  ~simulator(){}
  simulator(boost::property_tree::ptree &para_tree);

  int default_simulate(boost::property_tree::ptree &para_tree);
  int simulate(boost::property_tree::ptree &para_tree);
  int configMaterial(object *myobject,const double (&material_para)[2]); 
  int configForce(object *myobject);
};
#endif
