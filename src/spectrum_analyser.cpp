#include "spectrum_analyser.h"
#include "io.h"
#include "myvector.h"
#include "util.h"
using namespace std;
using namespace Eigen;

const static size_t vertex_index[8][3]=
  {
    {
      1,1,1
    },
    {
      1,0,1
    },
    {
      0,0,1
    },
    {
      0,1,1
    },
    {
      1,1,0
    },
    {
      1,0,0
    },
    {
      0,0,0
    },
    {
      0,1,0
    }
  };

spectrum_analyser::spectrum_analyser(boost::property_tree::ptree &para_tree)
{
  util myutil=util();
  io myio=io();
  string input_fine_object=para_tree.get<string>("input_fine_object.value");
  string input_fine_mat=para_tree.get<string>("input_fine_mat.value");
  string input_coarsen_object=para_tree.get<string>("input_coarsen_object.value");
  string input_coarsen_mat_new=para_tree.get<string>("input_coarsen_mat_new.value");

  string object_name=para_tree.get<string>("object_name.value");
  string out_dir=para_tree.get<string>("out_dir.value");
  int level=para_tree.get<int>("level.value"); //for edge: 1 to 2 ; 1 to 4 or 1 to 8

  size_t i,j,k;
  size_t a,b,c;
  size_t x,y,z;
  
  object* fine_obj=new object(input_fine_object,100,1,0,0,"co_rotated_linear"); 
  myio.getMatPara(fine_obj,input_fine_mat,1);

  object* coarsen_obj_new=new object(input_coarsen_object,100,1,0,0,"co_rotated_linear_with_stiffness_tensor"); 
  myio.getMatPara(coarsen_obj_new,input_coarsen_mat_new,0);
  
  size_t num_hex=coarsen_obj_new->myhexahedrons.size();
  size_t num_fine_in_coarsen=pow(level,3);
  
  for(i=0;i<0;++i) //coarsen hex index
    {
      MatrixXd spectrum=MatrixXd::Random(8*3,5);// fine, reileyCtoF, C, reileyC_oritoF, C_ori
      spectrum.fill(0);
      std::map<int,int > mpFromGlobalToLocal; mpFromGlobalToLocal.clear();
      std::map<int,int > mpFromLocalToGlobal; mpFromLocalToGlobal.clear();

      int num_vertex_local=0;
      for(a=0;a<2;a++)
	{
	  for(b=0;b<2;b++)
	    {
	      for(c=0;c<2;c++)
		{
		  int index_vertex_now=coarsen_obj_new->myhexahedrons[i].index_vertex[a][b][c];
		  mpFromGlobalToLocal[index_vertex_now]=num_vertex_local;
		  mpFromLocalToGlobal[num_vertex_local]=index_vertex_now;
		  num_vertex_local++;
		}
	    }
	}
      for(j=i*num_fine_in_coarsen;j<i*num_fine_in_coarsen+num_fine_in_coarsen;j++)
	{
	  for(a=0;a<2;a++)
	    {
	      for(b=0;b<2;b++)
		{
		  for(c=0;c<2;c++)
		    {
		      int index_vertex_now=fine_obj->myhexahedrons[j].index_vertex[a][b][c];
		      if(mpFromGlobalToLocal.find(index_vertex_now)!=mpFromGlobalToLocal.end())
			{
			  ;
			}
		      else
			{
			  mpFromGlobalToLocal[index_vertex_now]=num_vertex_local;
			  mpFromLocalToGlobal[num_vertex_local]=index_vertex_now;
			  num_vertex_local++;
			}
		    }
		}
	    }
	}
      printf("num_vertex_local: %d\n",num_vertex_local);

      MatrixXd Propagate_Matrix=MatrixXd::Random(num_vertex_local*3,8*3);
      Propagate_Matrix.fill(0);      
      calPropagate_Matrix(num_vertex_local,Propagate_Matrix,coarsen_obj_new,fine_obj,mpFromLocalToGlobal);
      cout<<"Propagate_Matrix:"<<endl;
      cout<<Propagate_Matrix<<endl;
      
      size_t num_vertex_fine=fine_obj->num_vertex; size_t num_vertex_coarsen=coarsen_obj_new->num_vertex;
      
      MatrixXd K_fine=MatrixXd::Random(num_vertex_fine*3,num_vertex_fine*3);  K_fine.fill(0);
      MatrixXd K_coarsen_ori=MatrixXd::Random(num_vertex_coarsen*3,num_vertex_coarsen*3); K_coarsen_ori.fill(0);
      MatrixXd K_coarsen_new=MatrixXd::Random(num_vertex_coarsen*3,num_vertex_coarsen*3); K_coarsen_new.fill(0);

      MatrixXd K_fine_local=MatrixXd::Random(num_vertex_local*3,num_vertex_local*3); K_fine_local.fill(0);
      MatrixXd K_coarsen_ori_local=MatrixXd::Random(8*3,8*3); K_coarsen_ori_local.fill(0);
      MatrixXd K_coarsen_new_local=MatrixXd::Random(8*3,8*3); K_coarsen_new_local.fill(0);

      MatrixXd Mass_fine_local=MatrixXd::Random(num_vertex_local*3,num_vertex_local*3); Mass_fine_local.fill(0);
      MatrixXd Mass_coarsen_ori_local=MatrixXd::Random(8*3,8*3); Mass_coarsen_ori_local.fill(0);
      MatrixXd Mass_coarsen_new_local=MatrixXd::Random(8*3,8*3); Mass_coarsen_new_local.fill(0);
      
      VectorXd J_fine(num_vertex_fine*3); J_fine.fill(0);
      VectorXd J_coarsen_ori(num_vertex_coarsen*3);  J_coarsen_ori.fill(0);
      VectorXd J_coarsen_new(num_vertex_coarsen*3);  J_coarsen_new.fill(0);

      MatrixXd y_new=MatrixXd::Random(8*3,8*3);
      MatrixXd y_ori=MatrixXd::Random(8*3,8*3);
      MatrixXd y_fine=MatrixXd::Random(num_vertex_local*3,num_vertex_local*3);

      MatrixXd X_mode_coarsen_new;
      //get K_coarsen_new_local      
      {
	coarsen_obj_new->myhexahedrons[i].calJacobianAndHessian(coarsen_obj_new->shapeFuncGrad,coarsen_obj_new->myvertexs,K_coarsen_new,J_coarsen_new,coarsen_obj_new->constitutive_model);
	// diff by hand 
	{
	  J_coarsen_new*=-1;
	}
	cout<<"K_coarsen_new:"<<endl;
	cout<<K_coarsen_new<<endl;
	// transform K_coarsen_new to local
	for(x=0;x<8;x++)
	  {
	    for(y=0;y<8;y++)
	      {
		int global_x=mpFromLocalToGlobal[x];
		int global_y=mpFromLocalToGlobal[y];
		for(a=0;a<3;a++)
		  {
		    for(b=0;b<3;b++)
		      {
			K_coarsen_new_local(3*x+a,3*y+b)=K_coarsen_new(3*global_x+a,3*global_y+b);
		      }
		  }
	      }
	  }
	cout<<"K_coarsen_new_local:"<<endl;
	cout<<K_coarsen_new_local<<endl;

	double vol_coarsen=coarsen_obj_new->myhexahedrons[i].avg_det_pX_peps*8;
	for(a=0;a<2;++a)
	  {
	    for(b=0;b<2;++b)
	      {
		for(c=0;c<2;++c)
		  {
		    x=mpFromGlobalToLocal[coarsen_obj_new->myhexahedrons[i].index_vertex[a][b][c]];
		    Mass_coarsen_new_local(3*x,3*x)+=(0.125*vol_coarsen*coarsen_obj_new->density);
		    Mass_coarsen_new_local(3*x+1,3*x+1)+=(0.125*vol_coarsen*coarsen_obj_new->density);
		    Mass_coarsen_new_local(3*x+2,3*x+2)+=(0.125*vol_coarsen*coarsen_obj_new->density);
		  }
	      }
	  }
	
	cout<<"Mass_coarsen_new_local's LLT's L matrix:"<<endl;
	MatrixXd L=Mass_coarsen_new_local.llt().matrixL();
	cout<<L<<endl;

	MatrixXd K_coarsen_new_use=L.inverse()*K_coarsen_new_local*(L.inverse().transpose());
       	SelfAdjointEigenSolver<MatrixXd> es1(K_coarsen_new_use);

        y_new=es1.eigenvectors();
	MatrixXd X_mode=L.transpose().inverse()*y_new;
	X_mode_coarsen_new=X_mode;
	spectrum.col(2)=es1.eigenvalues();


	coarsen_obj_new->myhexahedrons[i].ScalarPara=es1.eigenvalues();
	// check K
	
	for(j=0;j<24;j++)
	  {
	    stringstream jj; string strj; jj<<j; jj>>strj;
	    string out_local_mode_j=out_dir+"/"+object_name+"new."+strj+".vtk";
	    FILE*fp;
	    fp=fopen(out_local_mode_j.c_str(),"w");
	    fprintf(fp,"# vtk DataFile Version 2.0\n");
	    fprintf(fp,"tet\n");
	    fprintf(fp,"ASCII\n\n");
	    fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	    fprintf(fp,"POINTS 8 double\n");
	    for(a=0;a<8;a++)
	      {
       		fprintf(fp,"%lf %lf %lf\n",X_mode.col(j)[3*a]+coarsen_obj_new->myvertexs[mpFromLocalToGlobal[a]].location(0),X_mode.col(j)[3*a+1]+coarsen_obj_new->myvertexs[mpFromLocalToGlobal[a]].location(1),X_mode.col(j)[3*a+2]+coarsen_obj_new->myvertexs[mpFromLocalToGlobal[a]].location(2));
	      }
	    fprintf(fp,"CELLS 1 9\n");
	    for(a=0;a<1;++a)
	      {
		fprintf(fp,"8");
		for(b=0;b<8;++b)
		  {
		    fprintf(fp," %d",mpFromGlobalToLocal[coarsen_obj_new->myhexahedrons[i].index_vertex[vertex_index[b][0]][vertex_index[b][1]][vertex_index[b][2]]]);
		  }
		fprintf(fp,"\n");
	      }
	    fprintf(fp,"CELL_TYPES 1\n");
	    for(a=0;a<1;++a)
	      {
		fprintf(fp,"12\n");
	      }
	    fclose(fp);
	  }
      }
      
      //get K_fine_local
      
      {
	for(j=i*num_fine_in_coarsen;j<i*num_fine_in_coarsen+num_fine_in_coarsen;j++)
	  {
	    fine_obj->myhexahedrons[j].calJacobianAndHessian(fine_obj->shapeFuncGrad,fine_obj->myvertexs,K_fine,J_fine,fine_obj->constitutive_model);
	  }
	// diff by hand 
	{
	  J_fine*=-1;
	}
	for(x=0;x<num_vertex_local;x++)
	  {
	    for(y=0;y<num_vertex_local;y++)
	      {
		int global_x=mpFromLocalToGlobal[x];
		int global_y=mpFromLocalToGlobal[y];
		for(a=0;a<3;a++)
		  {
		    for(b=0;b<3;b++)
		      {
			K_fine_local(3*x+a,3*y+b)=K_fine(3*global_x+a,3*global_y+b);
		      }
		  }
	      }
	  }
	cout<<"K_fine_local:"<<endl;
	cout<<K_fine_local<<endl;

	for(j=i*num_fine_in_coarsen;j<i*num_fine_in_coarsen+num_fine_in_coarsen;j++)
	  {
	    double vol_fine=fine_obj->myhexahedrons[j].avg_det_pX_peps*8;
	    for(a=0;a<2;++a)
	      {
		for(b=0;b<2;++b)
		  {
		    for(c=0;c<2;++c)
		      {
			x=mpFromGlobalToLocal[fine_obj->myhexahedrons[j].index_vertex[a][b][c]];
			Mass_fine_local(3*x,3*x)+=(0.125*vol_fine*fine_obj->density);
			Mass_fine_local(3*x+1,3*x+1)+=(0.125*vol_fine*fine_obj->density);
			Mass_fine_local(3*x+2,3*x+2)+=(0.125*vol_fine*fine_obj->density);
		      }
		  }
	      }
	  }

	MatrixXd mass_diag=Mass_fine_local.diagonal();
	cout<<mass_diag<<endl;
	MatrixXd Res_mass=Propagate_Matrix.transpose()*mass_diag;
	cout<<"Res_mass:\n";
	cout<<Res_mass<<endl;
	
       	cout<<"Mass_fine_local's LLT's L matrix:"<<endl;
	MatrixXd L=Mass_fine_local.llt().matrixL();
	cout<<L<<endl;

	MatrixXd K_fine_use=L.inverse()*K_fine_local*(L.inverse().transpose());
	
	SelfAdjointEigenSolver<MatrixXd> es3(K_fine_use);
	y_fine=es3.eigenvectors();

	MatrixXd X_mode=L.transpose().inverse()*y_fine;

	spectrum.col(0)=es3.eigenvalues().segment(0,24);
	
	for(j=0;j<8*3;j++)
	  {
	    {
	      MatrixXd Pro_x_new(num_vertex_local*3,1); Pro_x_new=Propagate_Matrix*X_mode_coarsen_new.col(j);
	      MatrixXd help(1,1); help=(Pro_x_new.transpose()*Mass_fine_local*Pro_x_new);
	      MatrixXd help2(1,1); help2=(Pro_x_new.transpose()*K_fine_local*Pro_x_new)/help(0,0);
	      spectrum(j,1)=help2(0,0);
	    }
	  }

	
	for(j=0;j<8*3;j++)
	  {
	    stringstream jj; string strj; jj<<j; jj>>strj;
	    string out_local_mode_j=out_dir+"/"+object_name+"_fine."+strj+".vtk";
	    FILE*fp;
	    fp=fopen(out_local_mode_j.c_str(),"w");
	    fprintf(fp,"# vtk DataFile Version 2.0\n");
	    fprintf(fp,"tet\n");
	    fprintf(fp,"ASCII\n\n");
	    fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	    fprintf(fp,"POINTS %d double\n",num_vertex_local);

	    for(a=0;a<num_vertex_local;a++)
	      {
	       	fprintf(fp,"%lf %lf %lf\n",X_mode.col(j)[3*a]+fine_obj->myvertexs[mpFromLocalToGlobal[a]].location(0),X_mode.col(j)[3*a+1]+fine_obj->myvertexs[mpFromLocalToGlobal[a]].location(1),X_mode.col(j)[3*a+2]+fine_obj->myvertexs[mpFromLocalToGlobal[a]].location(2));
	      }
	    fprintf(fp,"CELLS %d %d\n",num_fine_in_coarsen,num_fine_in_coarsen*9);
	    for(a=0;a<num_fine_in_coarsen;++a)
	      {
		fprintf(fp,"8");
		for(b=0;b<8;++b)
		  {
		    fprintf(fp," %d",mpFromGlobalToLocal[fine_obj->myhexahedrons[i*num_fine_in_coarsen+a].index_vertex[vertex_index[b][0]][vertex_index[b][1]][vertex_index[b][2]]]);
		  }
		fprintf(fp,"\n");
	      }
	    fprintf(fp,"CELL_TYPES %d\n",num_fine_in_coarsen);
	    for(a=0;a<num_fine_in_coarsen;++a)
	      {
		fprintf(fp,"12\n");
	      }
	    fclose(fp);
	  }
      }
      cout<<"spectrum::\n";
      cout<<spectrum<<endl;
            
      {
	stringstream ss; string str; ss<<i; ss>>str;
	string out_local_spectrum=out_dir+"/"+object_name+"_spectrum."+str+".csv";
	FILE*fp;
	fp=fopen(out_local_spectrum.c_str(),"w");
	fprintf(fp,"fine,reileyCtoF,C,reileyC_oritoF,C_ori\n");
	for(a=0;a<24;a++)
	  {
	    for(b=0;b<5;b++)
	      {
		fprintf(fp,"%lf",spectrum(a,b));
		if(b!=5-1)
		  {
		    fprintf(fp,",");
		  }
	      }
	    fprintf(fp,"\n");
	  }
	fclose(fp);           
      }
    }

  MatrixXd spectrum_global=MatrixXd::Random(8*3,2);// fine,C
  spectrum_global.fill(0);

  MatrixXd spectrum_global_quick=MatrixXd::Random(36,2); //fine,C
  spectrum_global_quick.fill(0);
  //get K_fine_global
  {
    size_t num_vertex_fine=fine_obj->myvertexs.size();
    MatrixXd K_fine=MatrixXd::Random(num_vertex_fine*3,num_vertex_fine*3); K_fine.fill(0);
    VectorXd Mass_fine(num_vertex_fine*3);
    VectorXd J_fine(num_vertex_fine*3); J_fine.fill(0);
    for(i=0;i<fine_obj->myhexahedrons.size();i++)
      {
	fine_obj->myhexahedrons[i].calJacobianAndHessian(fine_obj->shapeFuncGrad,fine_obj->myvertexs,K_fine,J_fine,fine_obj->constitutive_model);
      }
    // diff by hand 
    {
      J_fine*=-1;
    }
    for(i=0;i<num_vertex_fine;i++)
      {
	Mass_fine(3*i+0)=Mass_fine(3*i+1)=Mass_fine(3*i+2)=fine_obj->myvertexs[i].mass;
      }
    size_t num_all_dof=num_vertex_fine*3;
    VectorXd L(num_all_dof);
    VectorXd L_inverse(num_all_dof);
    for(i=0;i<num_all_dof;++i)
      {
	L(i)=sqrt(Mass_fine(i));
	L_inverse(i)=1.0/L(i);
      }
    SparseMatrix<double > K_fine_use_spa(num_vertex_fine*3,num_vertex_fine*3);
    vector< Triplet<double > > tripletsFor_K_fine_use_spa;

    double EPS=1e-10;
    for(i=0;i<num_all_dof;++i)
      {
	for(j=0;j<num_all_dof;++j)
	  {
	    double val=K_fine(i,j)*L_inverse(i)*L_inverse(j);
	    if(fabs(val)>= EPS)
	      {
		tripletsFor_K_fine_use_spa.emplace_back(i,j,val);
	      }
	  }
      }
    K_fine_use_spa.setFromTriplets(tripletsFor_K_fine_use_spa.begin(),tripletsFor_K_fine_use_spa.end());
    K_fine_use_spa.makeCompressed();
    size_t eigen_num=36; VectorXd eigen_value(36);
    myutil.get_smallest_eigen_value(K_fine_use_spa,eigen_num,eigen_value);
    spectrum_global_quick.col(0)=eigen_value;
    /*
    SelfAdjointEigenSolver<MatrixXd> es3(K_fine_use);
    spectrum_global.col(0)=es3.eigenvalues().segment(0,24);
    MatrixXd y_fine=es3.eigenvectors();
    MatrixXd X_mode=L.transpose().inverse()*y_fine;
    */
    /*
    for(j=0;j<8*3;j++)
      {
	stringstream jj; string strj; jj<<j; jj>>strj;
	string out_global_mode_j=out_dir+"/"+object_name+"_fine_global."+strj+".vtk";
	FILE*fp;
	fp=fopen(out_global_mode_j.c_str(),"w");
	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"tet\n");
	fprintf(fp,"ASCII\n\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp,"POINTS %d double\n",fine_obj->myvertexs.size());

	double lambda=0.001;
	for(a=0;a<fine_obj->myvertexs.size();a++)
	  {
	    fprintf(fp,"%lf %lf %lf\n",X_mode.col(j)[3*a]*lambda+fine_obj->myvertexs[a].location(0),X_mode.col(j)[3*a+1]*lambda+fine_obj->myvertexs[a].location(1),X_mode.col(j)[3*a+2]*lambda+fine_obj->myvertexs[a].location(2));
	  }
	fprintf(fp,"CELLS %d %d\n",fine_obj->myhexahedrons.size(),fine_obj->myhexahedrons.size()*9);
	for(a=0;a<fine_obj->myhexahedrons.size();++a)
	  {
	    fprintf(fp,"8");
	    for(b=0;b<8;++b)
	      {
		fprintf(fp," %d",fine_obj->myhexahedrons[a].index_vertex[vertex_index[b][0]][vertex_index[b][1]][vertex_index[b][2]]);
	      }
	    fprintf(fp,"\n");
	  }
	fprintf(fp,"CELL_TYPES %d\n",fine_obj->myhexahedrons.size());
	for(a=0;a<fine_obj->myhexahedrons.size();++a)
	  {
	    fprintf(fp,"12\n");
	  }
	fclose(fp);
	}*/
  }

  //get K_coarsen_global
  {
    size_t num_vertex_coarsen=coarsen_obj_new->myvertexs.size();
    MatrixXd K_coarsen_new=MatrixXd::Random(num_vertex_coarsen*3,num_vertex_coarsen*3);  K_coarsen_new.fill(0);
    VectorXd Mass_coarsen_new(num_vertex_coarsen*3);
    VectorXd J_coarsen_new(num_vertex_coarsen*3); J_coarsen_new.fill(0);
    for(i=0;i<coarsen_obj_new->myhexahedrons.size();i++)
      {
	coarsen_obj_new->myhexahedrons[i].calJacobianAndHessian(coarsen_obj_new->shapeFuncGrad,coarsen_obj_new->myvertexs,K_coarsen_new,J_coarsen_new,coarsen_obj_new->constitutive_model);
      }
    // diff by hand 
    {
      J_coarsen_new*=-1;
    }
    for(i=0;i<num_vertex_coarsen;i++)
      {
	Mass_coarsen_new(3*i+0)=Mass_coarsen_new(3*i+1)=Mass_coarsen_new(3*i+2)=coarsen_obj_new->myvertexs[i].mass;
      }

    size_t num_all_dof=num_vertex_coarsen*3;
    VectorXd L(num_all_dof);
    VectorXd L_inverse(num_all_dof);
    for(i=0;i<num_all_dof;++i)
      {
	L(i)=sqrt(Mass_coarsen_new(i));
	L_inverse(i)=1.0/L(i);
      }
    SparseMatrix<double > K_coarsen_new_use_spa(num_vertex_coarsen*3,num_vertex_coarsen*3);
    vector< Triplet<double > > tripletsFor_K_coarsen_new_use_spa;

    
    double EPS=1e-10;

    for(i=0;i<num_all_dof;++i)
      {
	for(j=0;j<num_all_dof;++j)
	  {
	    double val=K_coarsen_new(i,j)*L_inverse(i)*L_inverse(j);
	    if(fabs(val)>= EPS)
	      {
		tripletsFor_K_coarsen_new_use_spa.emplace_back(i,j,val);
	      }
	  }
      }
    K_coarsen_new_use_spa.setFromTriplets(tripletsFor_K_coarsen_new_use_spa.begin(),tripletsFor_K_coarsen_new_use_spa.end());
    K_coarsen_new_use_spa.makeCompressed();

    size_t eigen_num=36; VectorXd eigen_value(36);
    myutil.get_smallest_eigen_value(K_coarsen_new_use_spa,eigen_num,eigen_value);
    
    spectrum_global_quick.col(1)=eigen_value;
    
    /*
    SelfAdjointEigenSolver<MatrixXd> es3(K_coarsen_new_use);
    spectrum_global.col(1)=es3.eigenvalues().segment(0,24);*/
  }

  /*
  {
    string out_global_spectrum=out_dir+"/"+object_name+"_global_spectrum.csv";
    FILE*fp;
    fp=fopen(out_global_spectrum.c_str(),"w");
    fprintf(fp,"fine_global,C_global\n");
    for(a=0;a<24;a++)
      {
	for(b=0;b<2;b++)
	  {
	    fprintf(fp,"%lf",spectrum_global(a,b));
	    if(b!=2-1)
	      {
		fprintf(fp,",");
	      }
	  }
	fprintf(fp,"\n");
      }
    fclose(fp);           
    }*/

  {
    string out_global_spectrum=out_dir+"/"+object_name+"_global_spectrum_quick.csv";
    FILE*fp;
    fp=fopen(out_global_spectrum.c_str(),"w");
    fprintf(fp,"fine_global,C_global\n");
    for(a=0;a<36;a++)
      {
	for(b=0;b<2;b++)
	  {
	    fprintf(fp,"%lf",spectrum_global_quick(a,b));
	    if(b!=2-1)
	      {
		fprintf(fp,",");
	      }
	  }
	fprintf(fp,"\n");
      }
    fclose(fp);           
  }
    

  /*
  //选一个 element的stiffness tensor 作为全部的
  for(i=0;i<num_hex;i++)
    {
      //  coarsen_obj_new->myhexahedrons[i].stiffness_tensor=coarsen_obj_new->myhexahedrons[0].stiffness_tensor;
    }
  string out_mat_para_coarsen=out_dir+"/"+object_name+".sub0.1.mat";
  myio.saveMatPara(coarsen_obj_new,out_mat_para_coarsen,0);
  
  for(i=0;i<24;i++)
    {
      stringstream ss; string str; ss<<i; ss>>str;
      string out_vtk=out_dir+"/"+object_name+"_coarsen_eigen_value."+str+".vtk";
      myio.saveAsVTKwithScalarPara(coarsen_obj_new,i,out_vtk);
      }*/
  delete coarsen_obj_new;
  delete fine_obj;
}

int spectrum_analyser::calPropagate_Matrix(int num_vertex_local,Eigen::MatrixXd &Pro_M,object* coarsen_obj,object* fine_obj,std::map<int,int > &mpFromLocalToGlobal)
{
  int i,j,a,b;
  MatrixXd pro=MatrixXd::Random(8,1);
  MatrixXd A=MatrixXd::Random(8,8);
  for(i=0;i<8;i++)
    {
      A(i,0)=1;
      A(i,1)=coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(0);
      A(i,2)=coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(1);
      A(i,3)=coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);
      A(i,4)=coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(0)*coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(1);
      A(i,5)=coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(0)*coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);
      A(i,6)=coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(1)*coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);
      A(i,7)=coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(0)*coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(1)*coarsen_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);
    }
  
  MatrixXd A_inverse=MatrixXd::Random(8,8); A_inverse=A.inverse();
  MatrixXd A_inverse_transpose=MatrixXd::Random(8,8); A_inverse_transpose=A_inverse.transpose();
  MatrixXd x_now=MatrixXd::Random(8,1);
  for(i=0;i<num_vertex_local;i++)
    {
      x_now(0,0)=1;
      x_now(1,0)=fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(0);
      x_now(2,0)=fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(1);
      x_now(3,0)=fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);
      x_now(4,0)=fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(0)*fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(1);
      x_now(5,0)=fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(0)*fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);
      x_now(6,0)=fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(1)*fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);
      x_now(7,0)=fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(0)*fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(1)*fine_obj->myvertexs[mpFromLocalToGlobal[i]].location(2);

      
      pro=A_inverse_transpose*x_now;
      for(j=0;j<8;j++)
	{
	  Pro_M(3*i,3*j)=Pro_M(3*i+1,3*j+1)=Pro_M(3*i+2,3*j+2)=pro(j,0);
	}
    }
  return 0;
}
