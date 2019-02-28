#ifndef _MATE_FACTORY_
#define _MATE_FACTORY_

#include "head.h"
#include "hexahedron.h"
#include "object.h"
class mate_factory
{
 public:
  double dmetric;
  ~mate_factory(){}
  mate_factory(){}
  mate_factory(boost::property_tree::ptree &para_tree);
  int mate_A(object* fine_obj,double y_homo,double p);
  int mate_B(object* fine_obj,double (&y)[2],double p,int axis);
  int mate_C(object* fine_obj,double (&y)[2],double p,int axis);
  int mate_D(object* fine_obj,double (&y)[2],double p);
  int mate_E(object* fine_obj,double (&y)[2],double p,int axis);
  int mate_F(object* fine_obj,double (&y)[2],double p);
  int mate_G(object* fine_obj,double (&y)[2],double p,double N,double C);
  // level set function is sin(N*pi*x)*cos(N*pi*y) + sin(N*pi*z)*cos(N*pi*x) + sin(N*pi*y)*cos(N*pi*z) - C
  // where negative is hard and positive is soft
  int cal_mate(hexahedron &hex_now,double y,double p);
};

#endif
