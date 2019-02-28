#include "hexahedron.h"
#include "autodiff.h"
#include "util.h"
using namespace std;
using namespace Eigen;

DECLARE_DIFFSCALAR_BASE();
hexahedron::hexahedron(const size_t (&index_vertex)[2][2][2])
{
  size_t i,j,k;
  for(i=0;i<2;++i)
    {
      for(j=0;j<2;++j)
	{
	  for(k=0;k<2;++k)
	    {
	      this->index_vertex[i][j][k]=index_vertex[i][j][k];
	    }
	}
    }
  stiffness_tensor.fill(0); //because material_para每个分量都会被占据,不会带来未赋值就使用的结果，但是stiffness tensor 在测试过程中可能只是被部分填充非零数
}

// precomputing the 1.inverse_pX_peps 2.det_pX_peps 3.avg_det_pX_peps
int hexahedron::prepare(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],std::vector< vertex > &myvertexs)
{
  size_t a,b,c;
  size_t i,j,k;
  MatrixXd F(3,3);
  MatrixXd shapeFuncGradNow(3,1);
  size_t row,col;
  size_t now;
  
  MatrixXd IM=MatrixXd::Identity(3,3);
  util myutil=util();
  size_t count_quadrature=0;
  for(a=0;a<2;++a)
    {
      for(b=0;b<2;++b)
	{
	  for(c=0;c<2;++c)
	    {	      
	      for(row=0;row<3;++row)
		{
		  for(col=0;col<3;++col)
		    {
		      double value_now=0;
		      for(i=0;i<2;++i)
			{
			  for(j=0;j<2;++j)
			    {
			      for(k=0;k<2;++k)
				{
				  value_now+=myvertexs[index_vertex[i][j][k]].location_original(row)*shapeFuncGrad[i][j][k][a][b][c][col];
	        		}
			    }
			}
		      F(row,col)=value_now;
		    }
		}
	      det_pX_peps[count_quadrature]=F.determinant();
	      inverse_pX_peps[count_quadrature]=F.inverse();	       
	      count_quadrature++;
	    }
	}
    }

  avg_det_pX_peps=0;
  for(i=0;i<8;++i)
    {
      avg_det_pX_peps+=det_pX_peps[i];
    }
  avg_det_pX_peps*=0.125;
  // printf("avg_det_pX_peps:%lf\n",avg_det_pX_peps);
  return 0;
}

bool hexahedron::checkInversion(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],std::vector< vertex > &myvertexs,Eigen::VectorXd &Jacobian,const std::string &model)
{
  size_t a,b,c;
  size_t i,j,k;
  MatrixXd F(3,3);
  MatrixXd E(3,3);
  MatrixXd P(3,3);
  MatrixXd R(3,3);
  MatrixXd RTF(3,3);
  MatrixXd shapeFuncGradNow(3,1);
  MatrixXd temAns(3,1); //temporary answer

  double EdotE;
  double EPS=1e-5; //EPS is local constant to judge zero 
  double energy=0;
  size_t row,col;
  size_t now;
  
  MatrixXd IM=MatrixXd::Identity(3,3);

  util myutil=util();

  size_t count_quadrature=0;
  for(a=0;a<2;++a)
    {
      for(b=0;b<2;++b)
	{
	  for(c=0;c<2;++c)
	    {	      
	      for(row=0;row<3;++row)
		{
		  for(col=0;col<3;++col)
		    {
		      double value_now=0;
		      for(i=0;i<2;++i)
			{
			  for(j=0;j<2;++j)
			    {
			      for(k=0;k<2;++k)
				{
				  value_now+=myvertexs[index_vertex[i][j][k]].location_maybe(row)*shapeFuncGrad[i][j][k][a][b][c][col];
	        		}
			    }
			}
		      F(row,col)=value_now;
		    }
		}
	      F=F*inverse_pX_peps[count_quadrature];

	     
	      if(model=="co_rotated_linear")
		{		  
		  //co_rotated_linear
		  myutil.extractRotation(F,R);
		  // JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
		  //  R=svd.matrixU()*(svd.matrixV().transpose());
		  RTF=R.transpose()*F;
		  P=2*material_para[a][b][c][0]*(F-R)+material_para[a][b][c][1]*(RTF.trace()-3)*R;
		  E=R.transpose()*F-IM;
		  EdotE=0;
		  for(row=0;row<3;++row)
		    {
		      for(col=0;col<3;++col)
			{
			  EdotE+=E(row,col)*E(row,col);
			}
		    }
		  energy+=(material_para[a][b][c][0]*EdotE+0.5*material_para[a][b][c][1]*pow(E.trace(),2))*det_pX_peps[count_quadrature];  
		}	      
	      for(i=0;i<2;++i) 
		{
		  for(j=0;j<2;++j)
		    {
		      for(k=0;k<2;++k)
			{
			  for(now=0;now<3;++now)
			    {
			      shapeFuncGradNow(now,0)=shapeFuncGrad[i][j][k][a][b][c][now];
			    }			    
			  temAns=P*(inverse_pX_peps[count_quadrature].transpose())*shapeFuncGradNow*det_pX_peps[count_quadrature];
			  for(now=0;now<3;++now)
			    {
			      Jacobian(index_vertex[i][j][k]*3+now)+=temAns(now,0);
			    }			  
			}
		    }
		}
	      if(F.determinant()<EPS)
		{
		  //	    return true; //for flip 
		}
	      count_quadrature++;
	    }
	}
    }
  this->energy_maybe=energy;
  //  printf("energy_maybe : %lf\n",energy_maybe);
  return false;
}
bool hexahedron::checkInversion(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],std::vector< vertex > &myvertexs,const std::string &model)
{
  size_t a,b,c;
  size_t i,j,k;
 
  MatrixXd F(3,3);
  MatrixXd R(3,3); //for corotated linear 
  MatrixXd E(3,3);
  MatrixXd E_array(6,1); //because E is symmetric, it can be compressed to 6*1 
  size_t row,col;

  double EdotE;
  double EPS=1e-5; //EPS is local constant to judge zero 

  double energy=0;
  
  MatrixXd IM=MatrixXd::Identity(3,3);
  util myutil=util();

  size_t count_quadrature=0;
  for(a=0;a<2;++a)
    {
      for(b=0;b<2;++b)
	{
	  for(c=0;c<2;++c)
	    {
	      for(row=0;row<3;++row)
		{
		  for(col=0;col<3;++col)
		    {
		      double value_now=0;
		      for(i=0;i<2;++i)
			{
			  for(j=0;j<2;++j)
			    {
			      for(k=0;k<2;++k)
				{
				  value_now+=myvertexs[index_vertex[i][j][k]].location_maybe(row)*shapeFuncGrad[i][j][k][a][b][c][col];
				}
			    }
			}
		      F(row,col)=value_now;
		    }
		}
	      F=F*inverse_pX_peps[count_quadrature]; //deformation gradient
	      	      
	      if(model=="stvk")
		{
		  //stvk model     
		  E=0.5*(F.transpose()*F-IM);   
		  EdotE=0;
		  for(row=0;row<3;++row)
		    {
		      for(col=0;col<3;++col)
			{
			  EdotE+=E(row,col)*E(row,col);
			}
		    }
		  energy+=(material_para[a][b][c][0]*EdotE+0.5*material_para[a][b][c][1]*pow(E.trace(),2))*det_pX_peps[count_quadrature];
		}
	      else if(model=="neo_hookean")
		{
		  //neo model
		  MatrixXd FTF(3,3);
		  FTF=F.transpose()*F;
		  double i1=FTF.trace();
		  double det=F.determinant();
		  double lg=log(det);
		  energy+=(material_para[a][b][c][0]*(i1-3)*0.5-material_para[a][b][c][0]*lg+material_para[a][b][c][1]*pow(lg,2)*0.5)*det_pX_peps[count_quadrature];
		}
	      else if(model=="co_rotated_linear")
		{
		  //co_rotated_linear
		  
		   myutil.extractRotation(F,R);
		   // JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
		   // R=svd.matrixU()*(svd.matrixV().adjoint());
		  E=R.transpose()*F-IM;
		  EdotE=0;
		  for(row=0;row<3;++row)
		    {
		      for(col=0;col<3;++col)
			{
			  EdotE+=E(row,col)*E(row,col);
			}
		    }
		  energy+=(material_para[a][b][c][0]*EdotE+0.5*material_para[a][b][c][1]*pow(E.trace(),2))*det_pX_peps[count_quadrature];  
		}
	      else if(model=="linear")
		{
		  //linear
		  E=0.5*(F.transpose()+F)-IM;
		  EdotE=0;
		  for(row=0;row<3;++row)
		    {
		      for(col=0;col<3;++col)
			{
			  EdotE+=E(row,col)*E(row,col);
			}
		    }
		  energy+=(material_para[a][b][c][0]*EdotE+0.5*material_para[a][b][c][1]*pow(E.trace(),2))*det_pX_peps[count_quadrature];
		}
	      else if(model=="linear_with_stiffness_tensor")
		{
		  // linear_with_stiffness_tensor
		  E=0.5*(F.transpose()+F)-IM;
		  E_array(0,0)=E(0,0); E_array(1,0)=E(1,1); E_array(2,0)=E(2,2); E_array(3,0)=2*E(1,2); E_array(4,0)=2*E(2,0); E_array(5,0)=2*E(0,1);
		  for(row=0;row<6;++row)
		    {
		      for(col=row;col<6;++col)
			{
			  if(row==col)
			    {
			      energy+=(stiffness_tensor(row,col)*E_array(row,0)*E_array(col,0)*0.5)*det_pX_peps[count_quadrature];
			    }
			  else if(row!=col)
			    {
			      energy+=(stiffness_tensor(row,col)*E_array(row,0)*E_array(col,0))*det_pX_peps[count_quadrature];
			    }
			}
		    }
		}
	      else if(model=="co_rotated_linear_with_stiffness_tensor")
		{
		  //co_rotated_linear_with_stiffness_tensor

		  myutil.extractRotation(F,R);
		  // JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
		  //  R=svd.matrixU()*(svd.matrixV().adjoint());
		  E=R.transpose()*F-IM;
		  E_array(0,0)=E(0,0); E_array(1,0)=E(1,1); E_array(2,0)=E(2,2); E_array(3,0)=2*E(1,2); E_array(4,0)=2*E(2,0); E_array(5,0)=2*E(0,1);
		  for(row=0;row<6;++row)
		    {
		      for(col=row;col<6;++col)
			{
			  if(row==col)
			    {
			      energy+=(stiffness_tensor(row,col)*E_array(row,0)*E_array(col,0)*0.5)*det_pX_peps[count_quadrature];
			    }
			  else if(row!=col)
			    {
			      energy+=(stiffness_tensor(row,col)*E_array(row,0)*E_array(col,0))*det_pX_peps[count_quadrature];
			    }
			}
		    }
		}
	      else if(model=="stvk_with_stiffness_tensor")
		{
		  //stvk_with_stiffness_tensor
		  E=0.5*(F.transpose()*F-IM);
		  E_array(0,0)=E(0,0); E_array(1,0)=E(1,1); E_array(2,0)=E(2,2); E_array(3,0)=2*E(1,2); E_array(4,0)=2*E(2,0); E_array(5,0)=2*E(0,1);
		  for(row=0;row<6;++row)
		    {
		      for(col=row;col<6;++col)
			{
			  if(row==col)
			    {
			      energy+=(stiffness_tensor(row,col)*E_array(row,0)*E_array(col,0)*0.5)*det_pX_peps[count_quadrature];
			    }
			  else if(row!=col)
			    {
			      energy+=(stiffness_tensor(row,col)*E_array(row,0)*E_array(col,0))*det_pX_peps[count_quadrature];
			    }
			}
		    }
		}     
	      if(F.determinant()<EPS)
		{
		  //    return true; // for flip
		}
	      count_quadrature++;
	    }
	}
    }
  this->energy_maybe=energy;
  //  printf("energy_maybe : %lf\n",energy_maybe);
  return false; 
}

int hexahedron::calJacobianAndHessian(const double (&shapeFuncGrad)[2][2][2][2][2][2][3],vector< vertex > &myvertexs,Eigen::MatrixXd &Hessian,Eigen::VectorXd &Jacobian,const std::string &model)
{
  size_t a,b,c;
  size_t i,j,k;
  size_t ii,jj,kk;
  MatrixXd F(3,3);
  MatrixXd E(3,3);
  MatrixXd E_array(6,1); //because E is symmetric, it can be compressed to 6*1
  MatrixXd P_array(6,1); // for linear elasticity, P is also symmetric,so P_array is just P for linear elasticity;but for corotational & stvk, it is not but this array still can express the P`s symmetric component;
  MatrixXd P_sym(3,3); //store the P`s symmetric component;
  MatrixXd dE_array(6,1); //because E is symmetric, dE is still symmetric;
  MatrixXd dP_array(6,1); //store the dP`s symmetric component;
  MatrixXd dP_sym(3,3);
  MatrixXd P(3,3);
  MatrixXd R(3,3);
  MatrixXd RTF(3,3);
  MatrixXd shapeFuncGradNow(3,1);
  MatrixXd temAns(3,1); //temporary answer
  MatrixXd temAnsForHessian(3,3);
  Matrix<double,3,3> dF[3];
  for(i=0;i<3;++i)
    {
      dF[i].fill(0);
    }
  MatrixXd dE(3,3);
  MatrixXd dP(3,3);
  
  size_t row,col;
  size_t now,cc;
  MatrixXd IM=MatrixXd::Identity(3,3);

  util myutil=util();
  /*
  {
    typedef DScalar2<double,VectorXd, MatrixXd> DScalar; // use DScalar2 for calculating gradient and hessian and use DScalar1 for calculating gradient

    DScalar F_d[3][3];
    DScalar F_d_ori[3][3];
    DScalar R_d[3][3];
    DScalar E_d[3][3];
    VectorXd x(24);
    size_t ct=0;
    for(a=0;a<2;++a)
      {
	for(b=0;b<2;++b)
	  {
	    for(c=0;c<2;++c)
	      {
		for(row=0;row<3;++row)
		  {
		    x(ct*3+row)=myvertexs[index_vertex[a][b][c]].location(row);
		  }
		ct++;
	      }
	  }
      }

    DiffScalarBase::setVariableCount(24);

    DScalar x_d[24];
    DScalar energy_d=DScalar(0);
    for(i=0;i<24;i++)
      {
	x_d[i]=DScalar(i,x(i));
      }

    size count_quadrature=0;
    for(a=0;a<2;++a)
      {
	for(b=0;b<2;++b)
	  {
	    for(c=0;c<2;++c)
	      {
		for(row=0;row<3;++row)
		  {
		    for(col=0;col<3;++col)
		      {
			DScalar value_d(0);
			ct=0;
			for(i=0;i<2;++i)
			  {
			    for(j=0;j<2;++j)
			      {
				for(k=0;k<2;++k)
				  {
				    value_d+=x_d[ct*3+row]*shapeFuncGrad[i][j][k][a][b][c][col];
				    ct++;
				  }
			      }
			  }
			F_d[row][col]=value_d;
		      }
		  }
		  for(row=0;row<3;++row)
		  {
		    for(col=0;col<3;++col)
		    {
		      F_d_ori[row][col]=F_d[row][col];
		    }
		  }
		for(row=0;row<3;++row)
		  {
		    for(col=0;col<3;++col)
		      {
		        F_d[row][col]=F_d_ori[row][0]*inverse_pX_peps[count_quadrature](0,col)+F_d_ori[row][1]*inverse_pX_peps[count_quadrature](1,col)+F_d_ori[row][2]*inverse_pX_peps[count_quadrature](2,col);
			F(row,col)=F_d[row][col].getValue();
		      }
		  }
		if(model=="co_rotated_linear")
		  {
		    //co_rotated_linear
		    JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
		    MatrixXd sv_m(3,3); sv_m=(svd.matrixU().transpose())*F*svd.matrixV();
		    for(row=0;row<3;++row)
		      {
			sv_m(row,row)=1/sv_m(row,row);
		      }
		    //  cout<<sv_m<<endl<<endl;

		    MatrixXd vsvt(3,3); vsvt=svd.matrixV()*sv_m*(svd.matrixV().transpose());

		    for(row=0;row<3;++row)
		      {
			for(col=0;col<3;++col)
			  {
			    R_d[row][col]=DScalar(0);
			    for(size_t xx=0;xx<3;++xx)
			      {
				R_d[row][col]+=F_d[row][xx]*vsvt(xx,col);
			      }
			  }
		      }
		    
		    R=svd.matrixU()*(svd.matrixV().transpose());
		    MatrixXd RT(3,3); RT=R.transpose();
		    DScalar RTF_d[3][3];
		    DScalar RFT_d[3][3];
		    for(row=0;row<3;++row)
		      {
			for(col=0;col<3;++col)
			  {
			    RTF_d[row][col]=DScalar(0);
			    RFT_d[row][col]=DScalar(0);
			    for(size_t xx=0;xx<3;++xx)
			      {
				RTF_d[row][col]+=RT(row,xx)*F_d[xx][col];
				RFT_d[row][col]+=R(row,xx)*F_d[col][xx];
			      }
			  }
		      }
		    for(row=0;row<3;++row)
		      {
			for(col=0;col<3;++col)
			  {
			    E_d[row][col]=DScalar(0);
			     //12年siggraph course 公式
			       for(size_t xx=0;xx<3;++xx)
			       {
			       //	E_d[row][col]+=R_d[xx][row]*F_d[xx][col];
			       E_d[row][col]+=RT(row,xx)*F_d[xx][col];
			       }
			       E_d[row][col]-=IM(row,col);
			       // E_d[row][col]=F_d[row][col]-R_d[row][col];

			    //02年vitual material 公式
			    //    E_d[row][col]=0.5*(RTF_d[row][col]+RTF_d[col][row])-IM(row,col);

			    //05年Efficient, physically plausible finite elements 公式
			    // E_d[row][col]=0.5*(F_d[row][col]+F_d[col][row]-R(row,col)-R(col,row));

			    //yy的公式
			    E_d[row][col]=0.5*(RTF_d[row][col]+RFT_d[row][col])-IM(row,col);
			  }
		      }
		    DScalar EdotE_d(0);
		    for(row=0;row<3;++row)
		      {
			for(col=0;col<3;++col)
			  {
			    EdotE_d+=E_d[row][col]*E_d[row][col];
			  }
		      }
		    DScalar E_trace_d=E_d[0][0]+E_d[1][1]+E_d[2][2];
		    energy_d+=(material_para[a][b][c][0]*EdotE_d+0.5*material_para[a][b][c][1]*pow(E_trace_d,2))*det_pX_peps[count_quadrature];	
		  }
		else if(model=="linear")
		  {
		    //linear

		    for(row=0;row<3;++row)
		      {
			for(col=0;col<3;++col)
			  {
			    E_d[row][col]=0.5*(F_d[row][col]+F_d[col][row])-IM(row,col);
			  }
		      }
		    DScalar EdotE_d(0);
		    for(row=0;row<3;++row)
		      {
			for(col=0;col<3;++col)
			  {
			    EdotE_d+=E_d[row][col]*E_d[row][col];
			  }
		      }
		    DScalar E_trace_d=E_d[0][0]+E_d[1][1]+E_d[2][2];
		    energy_d+=(material_para[a][b][c][0]*EdotE_d+0.5*material_para[a][b][c][1]*pow(E_trace_d,2))*det_pX_peps[count_quadrature];
		  }
		  count_quadrature++;
	      }
	  }
      }
      
    MatrixXd grad(24,1);
    MatrixXd hes(24,24);
    grad=energy_d.getGradient();
    hes=energy_d.getHessian();
    ct=0; size_t ct2;
    for(i=0;i<2;++i)
      {
	for(j=0;j<2;++j)
	  {
	    for(k=0;k<2;++k)
	      {
		for(row=0;row<3;++row)
		  {
		    Jacobian(index_vertex[i][j][k]*3+row)+=grad(ct*3+row);
		    ct2=0;
		    for(ii=0;ii<2;++ii)
		      {
			for(jj=0;jj<2;++jj)
			  {
			    for(kk=0;kk<2;++kk)
			      {
				for(col=0;col<3;col++)
				  {
				    Hessian(index_vertex[i][j][k]*3+row,index_vertex[ii][jj][kk]*3+col)+=hes(ct*3+row,ct2*3+col);
				  }
				ct2++;
			      }
			  }
		      }
		  }
		ct++;		
	      }
	  }
      }
  }
  */
  size_t count_quadrature=0;
  
  for(a=0;a<2;++a)
    {
      for(b=0;b<2;++b)
	{
	  for(c=0;c<2;++c)
	    {
	      for(row=0;row<3;++row)
		{
		  for(col=0;col<3;++col)
		    {
		      double value_now=0;
		      for(i=0;i<2;++i)
			{
			  for(j=0;j<2;++j)
			    {
			      for(k=0;k<2;++k)
				{
				  value_now+=myvertexs[index_vertex[i][j][k]].location(row)*shapeFuncGrad[i][j][k][a][b][c][col];
	        		}
			    }
			}
		      F(row,col)=value_now;
		    }
		}
	      F=F*inverse_pX_peps[count_quadrature];
	      
	      MatrixXd FTinverse(3,3);
	      MatrixXd Finverse(3,3);
	      double det=F.determinant();
	      // printf("det of F: %lf\n",det);
	      if(model=="stvk")
		{
		  //stvk model	      
		  E=0.5*(F.transpose()*F-IM);    
		  P=F*(2*material_para[a][b][c][0]*E+material_para[a][b][c][1]*E.trace()*IM);
		}
	      else if(model=="neo_hookean")
		{
		  //neo model		  
		  FTinverse=F.transpose().inverse();
		  Finverse=F.inverse();	  
		  P=material_para[a][b][c][0]*(F-FTinverse)+material_para[a][b][c][1]*log(det)*FTinverse;
		}
	      else if(model=="co_rotated_linear")
		{		  
		  //co_rotated_linear
		  myutil.extractRotation(F,R);
		  //  JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
		  // R=svd.matrixU()*(svd.matrixV().transpose());
		  RTF=R.transpose()*F;
		  P=2*material_para[a][b][c][0]*(F-R)+material_para[a][b][c][1]*(RTF.trace()-3)*R;
		}
	      else if(model=="linear")
		{
		  //linear
		  P=material_para[a][b][c][0]*(F+F.transpose()-2*IM)+material_para[a][b][c][1]*(F.trace()-3)*IM;
		}
	      else if(model=="linear_with_stiffness_tensor")
		{
		  // linear_with_stiffness_tensor
		  E=0.5*(F.transpose()+F)-IM;
		  E_array(0,0)=E(0,0); E_array(1,0)=E(1,1); E_array(2,0)=E(2,2); E_array(3,0)=2*E(1,2); E_array(4,0)=2*E(2,0); E_array(5,0)=2*E(0,1);
		  for(row=0;row<6;++row)
		    {
		      double sum=0;
		      for(col=0;col<6;++col)
			{
			  sum+=stiffness_tensor(row,col)*E_array(col,0);
			}
		      P_array(row,0)=sum;
		    }
		  P_sym(0,0)=P_array(0,0); P_sym(1,1)=P_array(1,0); P_sym(2,2)=P_array(2,0); P_sym(2,1)=P_sym(1,2)=P_array(3,0); P_sym(0,2)=P_sym(2,0)=P_array(4,0); P_sym(1,0)=P_sym(0,1)=P_array(5,0);
		  P=P_sym;
		}
	      else if(model=="co_rotated_linear_with_stiffness_tensor")
		{
		  //co_rotated_linear_with_stiffness_tensor
		   myutil.extractRotation(F,R);
		   // JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
		   // R=svd.matrixU()*(svd.matrixV().adjoint());
		  E=R.transpose()*F-IM;
		  E_array(0,0)=E(0,0); E_array(1,0)=E(1,1); E_array(2,0)=E(2,2); E_array(3,0)=2*E(1,2); E_array(4,0)=2*E(2,0); E_array(5,0)=2*E(0,1);
		  for(row=0;row<6;++row)
		    {
		      double sum=0;
		      for(col=0;col<6;++col)
			{
			  sum+=stiffness_tensor(row,col)*E_array(col,0);
			}
		      P_array(row,0)=sum;
		    }
		  P_sym(0,0)=P_array(0,0); P_sym(1,1)=P_array(1,0); P_sym(2,2)=P_array(2,0); P_sym(2,1)=P_sym(1,2)=P_array(3,0); P_sym(0,2)=P_sym(2,0)=P_array(4,0); P_sym(1,0)=P_sym(0,1)=P_array(5,0);
		  P=R*P_sym;
		}
	      else if(model=="stvk_with_stiffness_tensor")
		{
		  //stvk_with_stiffness_tensor
		  E=0.5*(F.transpose()*F-IM);
		  E_array(0,0)=E(0,0); E_array(1,0)=E(1,1); E_array(2,0)=E(2,2); E_array(3,0)=2*E(1,2); E_array(4,0)=2*E(2,0); E_array(5,0)=2*E(0,1);
		  for(row=0;row<6;++row)
		    {
		      double sum=0;
		      for(col=0;col<6;++col)
			{
			  sum+=stiffness_tensor(row,col)*E_array(col,0);
			}
		      P_array(row,0)=sum;
		    }
		  P_sym(0,0)=P_array(0,0); P_sym(1,1)=P_array(1,0); P_sym(2,2)=P_array(2,0); P_sym(2,1)=P_sym(1,2)=P_array(3,0); P_sym(0,2)=P_sym(2,0)=P_array(4,0); P_sym(1,0)=P_sym(0,1)=P_array(5,0);
		  P=F*P_sym;
		}
	      
	      for(i=0;i<2;++i) 
		{
		  for(j=0;j<2;++j)
		    {
		      for(k=0;k<2;++k)
			{
			  for(now=0;now<3;++now)
			    {
			      shapeFuncGradNow(now,0)=shapeFuncGrad[i][j][k][a][b][c][now];
			    }			    
			  temAns=P*(inverse_pX_peps[count_quadrature].transpose())*shapeFuncGradNow*det_pX_peps[count_quadrature];
			  for(now=0;now<3;++now)
			    {
			      Jacobian(index_vertex[i][j][k]*3+now)+=temAns(now,0);
			    }			  
			  for(ii=0;ii<2;++ii)
			    {
			      for(jj=0;jj<2;++jj)
				{
				  for(kk=0;kk<2;++kk)
				    {
				      for(now=0;now<3;++now)
					{
					  for(cc=0;cc<3;++cc)
					    {					      
						      dF[now](now,cc)=shapeFuncGrad[ii][jj][kk][a][b][c][0]*inverse_pX_peps[count_quadrature](0,cc)+shapeFuncGrad[ii][jj][kk][a][b][c][1]*inverse_pX_peps[count_quadrature](1,cc)+shapeFuncGrad[ii][jj][kk][a][b][c][2]*inverse_pX_peps[count_quadrature](2,cc);
					    }
					  if(model=="stvk")
					    {
					      //stvk model	  
					      dE=(dF[now].transpose()*F+F.transpose()*dF[now])*0.5;
					      dP=dF[now]*(2*material_para[a][b][c][0]*E+material_para[a][b][c][1]*E.trace()*IM)+F*(2*material_para[a][b][c][0]*dE+material_para[a][b][c][1]*dE.trace()*IM);
					  
					    }
					  else if(model=="neo_hookean")
					    {
					      //neo model
					      MatrixXd FinversedF(3,3);
					      FinversedF=Finverse*dF[now];
					      dP=material_para[a][b][c][0]*dF[now]+(material_para[a][b][c][0]-material_para[a][b][c][1]*log(det))*FTinverse*(dF[now].transpose())*FTinverse+material_para[a][b][c][1]*FinversedF.trace()*FTinverse;
					    }
					  else if(model=="co_rotated_linear")
					    {
					      //co_rotated_linear
					      dP=material_para[a][b][c][0]*(dF[now]+dF[now].transpose())+material_para[a][b][c][1]*dF[now].trace()*IM;	       
					    }
					  else if(model=="linear")
					    {
					      //linear
					      dP=material_para[a][b][c][0]*(dF[now]+dF[now].transpose())+material_para[a][b][c][1]*dF[now].trace()*IM; 
					    }
					  else if(model=="linear_with_stiffness_tensor")
					    {
					      //linear_with_stiffness_tensor
					      dE=(dF[now].transpose()+dF[now])*0.5;
					      dE_array(0,0)=dE(0,0); dE_array(1,0)=dE(1,1); dE_array(2,0)=dE(2,2); dE_array(3,0)=2*dE(1,2); dE_array(4,0)=2*dE(2,0); dE_array(5,0)=2*dE(0,1);
					      for(row=0;row<6;++row)
						{
						  double sum=0;
						  for(col=0;col<6;++col)
						    {
						      sum+=stiffness_tensor(row,col)*dE_array(col,0);
						    }
						  dP_array(row,0)=sum;
						}
					      dP_sym(0,0)=dP_array(0,0); dP_sym(1,1)=dP_array(1,0); dP_sym(2,2)=dP_array(2,0); dP_sym(2,1)=dP_sym(1,2)=dP_array(3,0); dP_sym(0,2)=dP_sym(2,0)=dP_array(4,0); dP_sym(1,0)=dP_sym(0,1)=dP_array(5,0);
					      dP=dP_sym;
					    }
					  else if(model=="co_rotated_linear_with_stiffness_tensor")
					    {
					      //co_rotated_linear_with_stiffness_tensor
					      dE=(dF[now].transpose()+dF[now])*0.5;
					      dE_array(0,0)=dE(0,0); dE_array(1,0)=dE(1,1); dE_array(2,0)=dE(2,2); dE_array(3,0)=2*dE(1,2); dE_array(4,0)=2*dE(2,0); dE_array(5,0)=2*dE(0,1);
					      for(row=0;row<6;++row)
						{
						  double sum=0;
						  for(col=0;col<6;++col)
						    {
						      sum+=stiffness_tensor(row,col)*dE_array(col,0);
						    }
						  dP_array(row,0)=sum;
						}
					      dP_sym(0,0)=dP_array(0,0); dP_sym(1,1)=dP_array(1,0); dP_sym(2,2)=dP_array(2,0); dP_sym(2,1)=dP_sym(1,2)=dP_array(3,0); dP_sym(0,2)=dP_sym(2,0)=dP_array(4,0); dP_sym(1,0)=dP_sym(0,1)=dP_array(5,0);
					      dP=dP_sym;
					    }
					  else if(model=="stvk_with_stiffness_tensor")
					    {
					      //stvk_with_stiffness_tensor
					      dE=(dF[now].transpose()*F+F.transpose()*dF[now])*0.5;
					      dE_array(0,0)=dE(0,0); dE_array(1,0)=dE(1,1); dE_array(2,0)=dE(2,2); dE_array(3,0)=2*dE(1,2); dE_array(4,0)=2*dE(2,0); dE_array(5,0)=2*dE(0,1);
					      for(row=0;row<6;++row)
						{
						  double sum=0;
						  for(col=0;col<6;++col)
						    {
						      sum+=stiffness_tensor(row,col)*dE_array(col,0);
						    }
						  dP_array(row,0)=sum;
						}
					      dP_sym(0,0)=dP_array(0,0); dP_sym(1,1)=dP_array(1,0); dP_sym(2,2)=dP_array(2,0); dP_sym(2,1)=dP_sym(1,2)=dP_array(3,0); dP_sym(0,2)=dP_sym(2,0)=dP_array(4,0); dP_sym(1,0)=dP_sym(0,1)=dP_array(5,0);
					      dP=dF[now]*P_sym+F*dP_sym;
					    }
					  temAns=dP*(inverse_pX_peps[count_quadrature].transpose())*shapeFuncGradNow*det_pX_peps[count_quadrature];
					  for(cc=0;cc<3;++cc)
					    {
					      temAnsForHessian(cc,now)=temAns(cc,0);
					    }				
					}
				      if(model=="co_rotated_linear"||model=="co_rotated_linear_with_stiffness_tensor")
					{
					  temAnsForHessian=R*temAnsForHessian*R.transpose();
					}				   
				      for(cc=0;cc<3;++cc)
					{
					  for(now=0;now<3;++now)
					    {
					       Hessian(index_vertex[i][j][k]*3+cc,index_vertex[ii][jj][kk]*3+now)+=temAnsForHessian(cc,now);
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	      count_quadrature++;
	    }
	}
    }
  return 0;
}
