#ifndef _OBJECT_
#define _OBJECT_

#include "head.h"
#include "myvector.h"
#include "vertex.h"
#include "hexahedron.h"
#include "face.h"

class object
{
 public:
  
  double time_all;
  double norm_Jacobian_cal;

  std::vector<std::pair<double,double > > time_norm_pair;

  std::vector<vertex > myvertexs; //need to save in vtk
  std::vector<hexahedron > myhexahedrons; //need to save in vtk
  std::vector<face > myfaces;
  size_t num_all_dof; // dof for all
  size_t num_cal_dof; // dof for calculating
  
  size_t num_vertex;
  size_t num_faces;
  size_t num_hexahedrons;
  
  size_t num_fixed;
  
  double length,width,height; 
  size_t*** index_for_vertex;
  size_t lc,wc,hc;

  size_t* mapIndexToLocInMartix;
  size_t* mapLocInMatrixToIndex;

  double shapeFuncGrad[2][2][2][2][2][2][3];
  
  bool converge;

  double density;
  //dt decides whether it is static or dynamic simulation 
  double dt;
  
  int line_search;
  double weight_line_search;
  size_t iteration_num;
  size_t max_iteration_num;
  std::string constitutive_model;

  std::map< std::pair<int, std::pair<int,int > > ,int> mpFromIndexVertexToIndexFace;//一个面由索引值最小的三个顶点的索引建立mapping
  myvector center_loc; //1.fix the zeroth and first moment 2.对于拓扑同胚于disk的obj方便判断外表面法线的朝向
  
  object();
  ~object();
  object(std::string input_dir,double dt,double density,int line_search,double weight_line_search,std::string constitutive_model);

  // for this mesh generation,we need the dmetric
  object(double dmetric,size_t l_size,size_t w_size,size_t h_size,double dt,double density,int line_search,double weight_line_search,std::string constitutive_model);

  int hex_prepare();
  int prepare();
  int checkFixedOrFree();
  int calShapeFuncGrad();
  int init_Energy_now_ForHex();
  int calMassForVertex();
   
  int dynamicSimulator(); // if d1dt equals zero,it is a static simulator also.
  int harmonic_def_static();
  
  int calJacobianAndHessianForHex(Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian);
  bool checkInversion(Eigen::VectorXd &dx);
  bool checkInversion(Eigen::VectorXd &dx,Eigen::VectorXd &Jacobian);
  bool checkInversion_for_harmonic_def(Eigen::VectorXd &dx);
  double calEnergyDif();
  double calElasticEnergy();
  int solve(Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian);
  int solve_harmonic_def_static(Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian);
};
#endif
