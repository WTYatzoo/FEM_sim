#include "coarsener.h"
#include "vertex.h"
#include "hexahedron.h"
#include "io.h"
using namespace std;
using namespace Eigen;
coarsener::coarsener(boost::property_tree::ptree &para_tree)
{
  coarsen_eq5(para_tree);
  coarsen_eq8(para_tree);
}
int coarsener::coarsen_eq5(boost::property_tree::ptree &para_tree)
{
  io myio=io();
  string input_fine_object=para_tree.get<string>("input_fine_object.value");
  string input_fine_mat=para_tree.get<string>("input_fine_mat.value");
  string input_coarsen_object=para_tree.get<string>("input_coarsen_object.value");
  string input_harmonic_def=para_tree.get<string>("input_harmonic_def.value"); // harmonic displacement 的下标从1开始，不是从0开始
  string object_name=para_tree.get<string>("object_name.value");
  string out_dir=para_tree.get<string>("out_dir.value");

  int level=para_tree.get<int>("level.value"); //for edge: 1 to 2 ; 1 to 4 or 1 to 8
  
  size_t i,j;
  size_t a,b;

  vector<vertex > myvertexs[6]; 
  vector<hexahedron > myhexahedrons[6];
  for(i=0;i<6;i++)
    {
      while(!myvertexs[i].empty())
	{
	  myvertexs[i].pop_back();
	}
      while(!myhexahedrons[i].empty())
	{
	  myhexahedrons[i].pop_back();
	}
    }
  size_t ct=0;
  for(i=0;i<3;++i)
    {
      for(j=i;j<3;++j)
	{
	  stringstream ss_i; string str_i; ss_i<<i+1; ss_i>>str_i;
	  stringstream ss_j; string str_j; ss_j<<j+1; ss_j>>str_j;
	  string path=input_harmonic_def+"_"+str_i+"_"+str_j+".vtk";
	  myio.getVertexAndHex(myvertexs[ct],myhexahedrons[ct],path);
       	  ct++;
	}
    }
  object* fine_obj=new object();
  myio.getVertexAndHex(fine_obj->myvertexs,fine_obj->myhexahedrons,input_fine_object);
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  myio.getMatPara(fine_obj,input_fine_mat,1);

  //转成 stiffness_tensor 形式
  for(i=0;i<fine_obj->num_hexahedrons;++i)
    {
      fine_obj->myhexahedrons[i].stiffness_tensor(0,0)=fine_obj->myhexahedrons[i].stiffness_tensor(1,1)=fine_obj->myhexahedrons[i].stiffness_tensor(2,2)=fine_obj->myhexahedrons[i].material_para[0][0][0][0]*2+fine_obj->myhexahedrons[i].material_para[0][0][0][1];
      fine_obj->myhexahedrons[i].stiffness_tensor(0,1)=fine_obj->myhexahedrons[i].stiffness_tensor(0,2)=fine_obj->myhexahedrons[i].stiffness_tensor(1,2)=fine_obj->myhexahedrons[i].stiffness_tensor(1,0)=fine_obj->myhexahedrons[i].stiffness_tensor(2,0)=fine_obj->myhexahedrons[i].stiffness_tensor(2,1)=fine_obj->myhexahedrons[i].material_para[0][0][0][1];
      fine_obj->myhexahedrons[i].stiffness_tensor(3,3)=fine_obj->myhexahedrons[i].stiffness_tensor(4,4)=fine_obj->myhexahedrons[i].stiffness_tensor(5,5)=fine_obj->myhexahedrons[i].material_para[0][0][0][0];
    }
  
  object* coarsen_obj=new object();
  myio.getVertexAndHex(coarsen_obj->myvertexs,coarsen_obj->myhexahedrons,input_coarsen_object);
  coarsen_obj->num_hexahedrons=coarsen_obj->myhexahedrons.size();
  coarsen_obj->num_vertex=coarsen_obj->myvertexs.size();

  double dmetric_coarsen=myio.getDmetric(input_coarsen_object);
  double dmetric_fine=dmetric_coarsen/(double)level;
  printf("dmetric_coarsen:: %lf\ndmetric_fine:: %lf\n",dmetric_coarsen,dmetric_fine);
  
  MatrixXd A=MatrixXd::Random(21,21);
  MatrixXd ATA=MatrixXd::Random(21,21);
  VectorXd B(21); MatrixXd B_m(21,1); MatrixXd ATB_m(21,1); VectorXd ATB(21);
  VectorXd C(21);
  size_t num_hex=coarsen_obj->myhexahedrons.size();

  size_t num_fine_in_coarsen=pow(level,3);
  printf("coarsen num_hex:: %d\n",num_hex);
  for(i=0;i<num_hex;++i) //coarsen hex index
    {
      size_t ct_21=0; //处理的当前行号
      A.fill(0); B.fill(0);
      for(a=0;a<6;++a) // harmonic dis index
	{
	  for(b=a;b<6;++b) // harmonic dis index
	    {  
	      for(j=num_fine_in_coarsen*i;j<num_fine_in_coarsen*i+num_fine_in_coarsen;++j) // fine hex index
		{
		  get_B(dmetric_fine,B,ct_21,myvertexs[a],myvertexs[b],fine_obj->myhexahedrons[j],fine_obj->shapeFuncGrad);
		}
	      get_A(dmetric_coarsen,A,ct_21,myvertexs[a],myvertexs[b],coarsen_obj->myhexahedrons[i],coarsen_obj->shapeFuncGrad);   
	      ct_21++;
	    }
	}
      A*=(pow(dmetric_coarsen,3)*0.125);
      ATA=A.transpose()*A;
      for(a=0;a<21;++a)
	{
	  B_m(a,0)=B(a);
	}
      ATB_m=A.transpose()*B_m;
      for(a=0;a<21;++a)
	{
	  ATB(a)=ATB_m(a,0);
	}

      SelfAdjointEigenSolver<MatrixXd> es(ATA);
      printf("condition number: %lf\n",fabs(es.eigenvalues()[20]/es.eigenvalues()[0]));
      
      C=ATA.llt().solve(ATB);
      //  C=A.lu().solve(B); //the same answer as the previous line`s answer


      
      
      size_t ct_c=0;
      for(a=0;a<6;++a)
	{
	  for(b=a;b<6;++b)
	    {
	      coarsen_obj->myhexahedrons[i].stiffness_tensor(a,b)=coarsen_obj->myhexahedrons[i].stiffness_tensor(b,a)=C(ct_c);
	      ct_c++;
	    }
	}
      cout<<coarsen_obj->myhexahedrons[i].stiffness_tensor<<endl;
    }
  i=0;
  stringstream ss; string str; ss<<i; ss>>str;
  string out_mat_para_coarsen=out_dir+"/"+object_name+".sub"+str+".mat";
  
  myio.saveMatPara(coarsen_obj,out_mat_para_coarsen,0);

  delete coarsen_obj;
  delete fine_obj;
}


int coarsener::coarsen_eq8(boost::property_tree::ptree &para_tree)
{
  io myio=io();
  string input_fine_object=para_tree.get<string>("input_fine_object.value");
  string input_fine_mat=para_tree.get<string>("input_fine_mat.value");
  string input_coarsen_object=para_tree.get<string>("input_coarsen_object.value");
  string input_harmonic_def=para_tree.get<string>("input_harmonic_def.value"); // harmonic displacement 的下标从1开始，不是从0开始
  string object_name=para_tree.get<string>("object_name.value");
  string out_dir=para_tree.get<string>("out_dir.value");

  int level=para_tree.get<int>("level.value"); //for edge: 1 to 2 ; 1 to 4 or 1 to 8
  
  size_t i,j,k;
  size_t a,b,c;

  vector<vertex > myvertexs[6]; 
  vector<hexahedron > myhexahedrons[6];
  for(i=0;i<6;i++)
    {
      while(!myvertexs[i].empty())
	{
	  myvertexs[i].pop_back();
	}
      while(!myhexahedrons[i].empty())
	{
	  myhexahedrons[i].pop_back();
	}
    }
  size_t ct=0;
  for(i=0;i<3;++i)
    {
      for(j=i;j<3;++j)
	{
	  stringstream ss_i; string str_i; ss_i<<i+1; ss_i>>str_i;
	  stringstream ss_j; string str_j; ss_j<<j+1; ss_j>>str_j;
	  string path=input_harmonic_def+"_"+str_i+"_"+str_j+".vtk";
	  myio.getVertexAndHex(myvertexs[ct],myhexahedrons[ct],path);
       	  ct++;
	}
    }
  object* fine_obj=new object();
  myio.getVertexAndHex(fine_obj->myvertexs,fine_obj->myhexahedrons,input_fine_object);
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  myio.getMatPara(fine_obj,input_fine_mat,1);

  //转成 stiffness_tensor 形式
  for(i=0;i<fine_obj->num_hexahedrons;++i)
    {
      fine_obj->myhexahedrons[i].stiffness_tensor(0,0)=fine_obj->myhexahedrons[i].stiffness_tensor(1,1)=fine_obj->myhexahedrons[i].stiffness_tensor(2,2)=fine_obj->myhexahedrons[i].material_para[0][0][0][0]*2+fine_obj->myhexahedrons[i].material_para[0][0][0][1];
      fine_obj->myhexahedrons[i].stiffness_tensor(0,1)=fine_obj->myhexahedrons[i].stiffness_tensor(0,2)=fine_obj->myhexahedrons[i].stiffness_tensor(1,2)=fine_obj->myhexahedrons[i].stiffness_tensor(1,0)=fine_obj->myhexahedrons[i].stiffness_tensor(2,0)=fine_obj->myhexahedrons[i].stiffness_tensor(2,1)=fine_obj->myhexahedrons[i].material_para[0][0][0][1];
      fine_obj->myhexahedrons[i].stiffness_tensor(3,3)=fine_obj->myhexahedrons[i].stiffness_tensor(4,4)=fine_obj->myhexahedrons[i].stiffness_tensor(5,5)=fine_obj->myhexahedrons[i].material_para[0][0][0][0];
    }
  
  object* coarsen_obj=new object();
  myio.getVertexAndHex(coarsen_obj->myvertexs,coarsen_obj->myhexahedrons,input_coarsen_object);
  coarsen_obj->num_hexahedrons=coarsen_obj->myhexahedrons.size();
  coarsen_obj->num_vertex=coarsen_obj->myvertexs.size();

  double dmetric_coarsen=myio.getDmetric(input_coarsen_object);
  double dmetric_fine=dmetric_coarsen/(double)level;
  printf("dmetric_coarsen:: %lf\ndmetric_fine:: %lf\n",dmetric_coarsen,dmetric_fine);

  size_t num_hex=coarsen_obj->myhexahedrons.size();
  size_t num_fine_in_coarsen=pow(level,3);
  printf("coarsen num_hex:: %d\n",num_hex);
  for(i=0;i<num_hex;++i) //coarsen hex index
    {
      MatrixXd C(6,6); C.fill(0);
      for(j=num_fine_in_coarsen*i;j<num_fine_in_coarsen*i+num_fine_in_coarsen;++j) // fine hex index
	{
	  MatrixXd GT_fine_now(6,6);
	  for(a=0;a<6;a++)
	    {
	      get_GT(dmetric_fine,GT_fine_now,a,myvertexs[a],fine_obj->myhexahedrons[j],fine_obj->shapeFuncGrad);
	    }
	  C+=(GT_fine_now*fine_obj->myhexahedrons[j].stiffness_tensor*(GT_fine_now.transpose()));
	}
      C*=(1.0/(double)num_fine_in_coarsen);
      MatrixXd GT_coarsen(6,6);
      for(a=0;a<6;a++)
	{
	  get_GT(dmetric_coarsen,GT_coarsen,a,myvertexs[a],coarsen_obj->myhexahedrons[i],coarsen_obj->shapeFuncGrad);   
	}
      MatrixXd GT_inverse_coarsen(6,6); GT_inverse_coarsen=GT_coarsen.inverse();
      MatrixXd C_final(6,6); C_final=GT_inverse_coarsen*C*(GT_inverse_coarsen.transpose());

      SelfAdjointEigenSolver<MatrixXd> es(C_final);
      printf("condition number equation8: %lf\n",fabs(es.eigenvalues()[5]/es.eigenvalues()[0]));
      
      //update coarsen_obj stiffness_tensor
      coarsen_obj->myhexahedrons[i].stiffness_tensor=C_final;
      cout<<coarsen_obj->myhexahedrons[i].stiffness_tensor<<endl;
    }
  
  string out_mat_para_coarsen=out_dir+"/"+object_name+".sub0.1.mat"; //numerical coarsening's equation 8's result
  myio.saveMatPara(coarsen_obj,out_mat_para_coarsen,0);
  delete coarsen_obj;
  delete fine_obj;
  
}

int coarsener::get_GT(const double &dmetric,Eigen::MatrixXd &GT,const size_t &which, std::vector< vertex > &myvertexs,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][2][2][2][3])
{
  size_t a,b,c;
  size_t i,j,k;
 
  MatrixXd F(3,3);
  MatrixXd E(3,3);
  MatrixXd E_array(1,6);
  
  size_t row,col;
  
  double help=2.0/dmetric;
  MatrixXd IM=MatrixXd::Identity(3,3);

  for(a=0;a<2;++a) //quadrature point
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
		      for(i=0;i<2;++i) //顶点
			{
			  for(j=0;j<2;++j)
			    {
			      for(k=0;k<2;++k)
				{
				  value_now+=myvertexs[hexahedron_here.index_vertex[i][j][k]].location(row)*shapeFuncGrad[i][j][k][a][b][c][col];
				}
			      
			    }
			}
		      F(row,col)=value_now;
		    }
		}
	      F*=help; //deformation gradient	        
	      E=0.5*(F.transpose()+F)-IM;
	      E_array(0,0)=E(0,0); E_array(1,0)=E(1,1); E_array(2,0)=E(2,2); E_array(3,0)=2*E(1,2); E_array(4,0)=2*E(2,0); E_array(5,0)=2*E(0,1);

	      // here is a bug 
	      GT.row(which)=E_array;
	    }	  
	}
    }
  return 0;
}

int coarsener::get_B(const double &dmetric,Eigen::VectorXd &B,const size_t &ct_21, std::vector< vertex > &myvertexs1, std::vector< vertex > &myvertexs2,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][2][2][2][3])
{
  size_t a,b,c;
  size_t i,j,k;
 
  MatrixXd F1(3,3);
  MatrixXd E1(3,3);
  MatrixXd E_array1(6,1);
  MatrixXd F2(3,3);
  MatrixXd E2(3,3);
  MatrixXd E_array2(6,1);
  
  size_t row,col;
  double energy=0;

  double help=2.0/dmetric;
  MatrixXd IM=MatrixXd::Identity(3,3);

  for(a=0;a<2;++a) //quadrature point
    {
      for(b=0;b<2;++b)
	{
	  for(c=0;c<2;++c)
	    {
	      for(row=0;row<3;++row)
		{
		  for(col=0;col<3;++col)
		    {
		      double value_now1=0;
		      double value_now2=0;
		      for(i=0;i<2;++i) //顶点
			{
			  for(j=0;j<2;++j)
			    {
			      for(k=0;k<2;++k)
				{
				  value_now1+=myvertexs1[hexahedron_here.index_vertex[i][j][k]].location(row)*shapeFuncGrad[i][j][k][a][b][c][col];
				  value_now2+=myvertexs2[hexahedron_here.index_vertex[i][j][k]].location(row)*shapeFuncGrad[i][j][k][a][b][c][col];
				}
			    }
			}
		      F1(row,col)=value_now1;
		      F2(row,col)=value_now2;
		    }
		}
	      F1*=help; //deformation gradient
	      F2*=help;
	        
	      E1=0.5*(F1.transpose()+F1)-IM;
	      E2=0.5*(F2.transpose()+F2)-IM;
	      E_array1(0,0)=E1(0,0); E_array1(1,0)=E1(1,1); E_array1(2,0)=E1(2,2); E_array1(3,0)=2*E1(1,2); E_array1(4,0)=2*E1(2,0); E_array1(5,0)=2*E1(0,1);
	      E_array2(0,0)=E2(0,0); E_array2(1,0)=E2(1,1); E_array2(2,0)=E2(2,2); E_array2(3,0)=2*E2(1,2); E_array2(4,0)=2*E2(2,0); E_array2(5,0)=2*E2(0,1);

	      for(row=0;row<6;++row)
		{
		  for(col=row;col<6;++col)
		    {
		      if(row==col)
			{
			  energy+=hexahedron_here.stiffness_tensor(row,col)*E_array1(row,0)*E_array2(col,0);
			}
		      else if(row!=col)
			{
			  energy+=hexahedron_here.stiffness_tensor(row,col)*(E_array1(row,0)*E_array2(col,0)+E_array2(row,0)*E_array1(col,0));
			}
		    }
		}	
	    }
	}
    }
  energy*=(pow(dmetric,3)*0.125);
  B(ct_21)+=energy;
  return 0;
}
int coarsener::get_A(const double &dmetric,Eigen::MatrixXd &A,const size_t &ct_21, std::vector< vertex > &myvertexs1, std::vector< vertex > &myvertexs2,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][2][2][2][3])
{
  size_t a,b,c;
  size_t i,j,k;
 
  MatrixXd F1(3,3);
  MatrixXd E1(3,3);
  MatrixXd E_array1(6,1);
  MatrixXd F2(3,3);
  MatrixXd E2(3,3);
  MatrixXd E_array2(6,1);
  
  size_t row,col;
  
  double help=2.0/dmetric;
  MatrixXd IM=MatrixXd::Identity(3,3);

  for(a=0;a<2;++a) //quadrature point
    {
      for(b=0;b<2;++b)
	{
	  for(c=0;c<2;++c)
	    {
	      for(row=0;row<3;++row)
		{
		  for(col=0;col<3;++col)
		    {
		      double value_now1=0;
		      double value_now2=0;
		      for(i=0;i<2;++i) //顶点
			{
			  for(j=0;j<2;++j)
			    {
			      for(k=0;k<2;++k)
				{
				  value_now1+=myvertexs1[hexahedron_here.index_vertex[i][j][k]].location(row)*shapeFuncGrad[i][j][k][a][b][c][col];
				  value_now2+=myvertexs2[hexahedron_here.index_vertex[i][j][k]].location(row)*shapeFuncGrad[i][j][k][a][b][c][col];
				}
			    }
			}
		      F1(row,col)=value_now1;
		      F2(row,col)=value_now2;
		    }
		}
	      F1*=help; //deformation gradient
	      F2*=help;
	      
       	      E1=0.5*(F1.transpose()+F1)-IM;
	      E2=0.5*(F2.transpose()+F2)-IM;
	      E_array1(0,0)=E1(0,0); E_array1(1,0)=E1(1,1); E_array1(2,0)=E1(2,2); E_array1(3,0)=2*E1(1,2); E_array1(4,0)=2*E1(2,0); E_array1(5,0)=2*E1(0,1);
	      E_array2(0,0)=E2(0,0); E_array2(1,0)=E2(1,1); E_array2(2,0)=E2(2,2); E_array2(3,0)=2*E2(1,2); E_array2(4,0)=2*E2(2,0); E_array2(5,0)=2*E2(0,1);
	      
	      size_t ct_col=0;
	      for(row=0;row<6;++row)
		{
		  for(col=row;col<6;++col)
		    {
		      if(row==col)
			{
			  A(ct_21,ct_col)+=(E_array1(row,0)*E_array2(col,0));
			}
		      else if(row!=col)
			{
			  A(ct_21,ct_col)+=(E_array1(row,0)*E_array2(col,0)+E_array2(row,0)*E_array1(col,0));
			}
		      ct_col++;
		    }
		}	
	    }
	}
    }
  
  return 0;
}
