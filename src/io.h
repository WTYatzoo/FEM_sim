#ifndef _IO_
#define _IO_
#include "head.h"
#include "vertex.h"
#include "hexahedron.h"
#include "object.h"
class io
{
 public:
  io(){}
  ~io(){}
  //name 不用引用，因为可能实参直接是"xxx" 
  int saveAsVTK(object* myobject,const std::string name);
  int saveAsVTKwithForce(object* myobject,const std::string name);
  int saveAsVTKwithPara(object* myobject,const std::string name);
  int saveAsVTKwithScalarPara(object* myobject,int which,const std::string name);
  int saveArrayAsVTK(std::vector<vertex> &myvertexs,std::vector<hexahedron> &myhexahedrons,const std::string name);
  int saveMatPara(object* myobject,const std::string name,int kind);
  int getMatPara(object* myobject,const std::string name,int kind);
  int getConstraintFromCsv(object* myobject,const std::string name);
  //for uniform lattice the following function can be used 
  double getDmetric(const std::string name);
  int getVertexAndHex(std::vector<vertex> &myvertexs,std::vector<hexahedron> &myhexahedrons,const std::string name);
  int getMatParaFromTXT(object* myobject,const std::string name);
  int saveMatParaAsTXT(object* myobject,const std::string name);
  int save_ERROR_TXT(double error_value,const std::string name);
  int saveTimeNormPair(object* myobject,const std::string name);
};
#endif
