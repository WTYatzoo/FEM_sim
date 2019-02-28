#ifndef _UTIL_
#define _UTIL_

#include "head.h"

class util
{
 public:
  util(){}
  ~util(){}
  int extractRotation(Eigen::MatrixXd &B,Eigen::MatrixXd &R);
  int get_smallest_eigen_value(Eigen::SparseMatrix<double> &K_use,size_t eigen_num,Eigen::VectorXd &eigenvalues);
  
};

#endif
