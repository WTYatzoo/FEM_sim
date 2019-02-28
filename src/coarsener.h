#ifndef _COARSENER_
#define _COARSENER_

#include "head.h"
#include "vertex.h"
#include "hexahedron.h"

class coarsener
{
 public:
  ~coarsener(){}
  coarsener(){}
  coarsener(boost::property_tree::ptree &para_tree);
  int coarsen_eq5(boost::property_tree::ptree &para_tree);
  int coarsen_eq8(boost::property_tree::ptree &para_tree);
  int get_B(const double &dmetric,Eigen::VectorXd &B,const size_t &ct_21, std::vector< vertex > &myvertexs1, std::vector< vertex > &myvertexs2,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][2][2][2][3]);
  int get_A(const double &dmetric,Eigen::MatrixXd &A,const size_t &ct_21, std::vector< vertex > &myvertexs1, std::vector< vertex > &myvertexs2,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][2][2][2][3]);
  int get_GT(const double &dmetric,Eigen::MatrixXd &GT,const size_t &which, std::vector< vertex > &myvertexs,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][2][2][2][3]);
};

#endif
