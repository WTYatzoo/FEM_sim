#include "util.h"
using namespace Eigen;
//using namespace Spectra;

int util::extractRotation(Eigen::MatrixXd &B,Eigen::MatrixXd &R_return)
{
  Matrix3d A; A=B;
  Quaterniond q=Quaterniond(A);
  q.normalize();
  for(size_t iter=0;iter<5;iter++)
    {      
      Matrix3d R=q.matrix();
      Vector3d omega=(R.col(0).cross(A.col(0))+R.col(1).cross(A.col(1))+R.col(2).cross(A.col(2)))*(1.0/fabs(R.col(0).dot(A.col(0))+R.col(1).dot(A.col(1))+R.col(2).dot(A.col(2)))+1e-9);
      double w=omega.norm();
      if(w<1e-9)
	{
	  break; 
	}
      q=Quaterniond(AngleAxisd(w,(1.0/w)*omega))*q;
      q.normalize();
    }
  R_return=q.matrix();

  return 0;
}

int util::get_smallest_eigen_value(Eigen::SparseMatrix<double> &K_use,size_t eigen_num,Eigen::VectorXd &eigenvalues)
{
  /*
  SparseSymMatProd<double> op(K_use);
  const int maxits = 20000;
  const double tolerance = 1e-10;
  const int ncv = min((int)K_use.cols(), 2*(int)eigen_num);
  Spectra::SymEigsSolver<double, SMALLEST_ALGE, SparseSymMatProd<double>> solver(&op, eigen_num, ncv);
  solver.init();
  int nconv = solver.compute(maxits, tolerance, SMALLEST_ALGE);
  eigenvalues = solver.eigenvalues();
  */
  return 0;
}
