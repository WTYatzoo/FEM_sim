#ifndef _HEXAHEDRON_
#define _HEXAHEDRON_

#include "head.h"
#include "vertex.h"
class hexahedron
{
 public:
  size_t index_vertex[2][2][2];
  double energy_now;
  double energy_maybe;
  Eigen::Matrix<double,24,1> ScalarPara;
  double material_para[2][2][2][2]; //material para at eight quadrature points
  Eigen::Matrix<double,6,6> stiffness_tensor; // stiffness tensor is rank-4 tensor but it can be reduced to rank-2 as 6*6 symmetric matrix

  //for general hex,we need patial X/patial epsilon on the 8 quadrature points . Often we use the inverse of the matrix,so we save the inverse.
  Eigen::Matrix<double,3,3> inverse_pX_peps[8];
  double det_pX_peps[8];
  double avg_det_pX_peps;
  
  //以上如果在element中是各项同性的可以简化用material_para表示材料属性，但是如果是拟合出来的不再是各项同性的将使用stiffness_tensor表示材料属性
  hexahedron(){}
  ~hexahedron(){}
  hexahedron(const size_t (&index_vertex)[2][2][2]);
  int prepare(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],std::vector< vertex > &myvertexs);
  int calJacobianAndHessian(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],std::vector< vertex > &myvertexs,Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian,const std::string &model);
  bool checkInversion(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],std::vector< vertex > &myvertexs,const std::string &model);
  bool checkInversion(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],std::vector< vertex > &myvertexs,Eigen::VectorXd &Jacobian,const std::string &model);
};
#endif
