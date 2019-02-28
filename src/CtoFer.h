#ifndef _CTOFER_
#define _CTOFER_
#include "head.h"
#include "vertex.h"
#include "hexahedron.h"
#include "object.h"
class CtoFer
{
 public:
  ~CtoFer(){}
  CtoFer(){}
  CtoFer(boost::property_tree::ptree &para_tree);
  int calH(Eigen::MatrixXd &col_H,std::vector<vertex > myvertexs,object *fine_obj,size_t index_vertex_here);
  int calHbc(Eigen::MatrixXd &col_H_bc, std::vector<vertex > &myvertexs,const std::vector<hexahedron> &myhexahedrons,size_t i,size_t num_fine_in_coarsen,Eigen::MatrixXd &bc_ori);
  int calHbc(Eigen::MatrixXd &col_H_bc,std::vector<vertex > &myvertexs,const hexahedron &hexahedron_here,Eigen::MatrixXd &bc_ori);
  int calBaryCenter(Eigen::MatrixXd &bc_ori,Eigen::MatrixXd &bc, std::vector<vertex > &myvertexs,const hexahedron &hexahedron_here);
  int calRandS(const double &dmetric,Eigen::MatrixXd &R,Eigen::MatrixXd &S, std::vector< vertex > &myvertexs,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][3]);
  int calShapeFuncGrad(double (&shapeFuncGrad)[2][2][2][3]);
};
#endif
