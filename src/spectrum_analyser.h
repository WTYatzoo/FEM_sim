#ifndef _SPECTRUM_ANALYSER_
#define _SPECTRUM_ANALYSER_

#include "head.h"
#include "object.h"
class spectrum_analyser
{
 public:
  spectrum_analyser(){}
  ~spectrum_analyser(){}
  spectrum_analyser(boost::property_tree::ptree &para_tree);
  int calPropagate_Matrix(int num_vertex_local,Eigen::MatrixXd &Pro_M,object* coarsen_obj,object* fine_obj,std::map<int,int > &mpFromLocalToGlobal);
};
#endif
