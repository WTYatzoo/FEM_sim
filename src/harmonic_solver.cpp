#include "harmonic_solver.h"
#include "io.h"
using namespace std;
using namespace Eigen;
harmonic_solver::harmonic_solver(boost::property_tree::ptree &para_tree)
{
  string kind=para_tree.get<string>("kind.value");
  scale=para_tree.get<double>("scale.value",100);
  if(kind=="default")
    {
      default_solve(para_tree);
    }
  else if(kind=="use_data_from_file")
    {
      solve(para_tree);
    }
}

int harmonic_solver::default_solve(boost::property_tree::ptree &para_tree)
{
  io myio=io();
  string out_dir=para_tree.get<string>("out_dir.value");
  double PoissonRatio=para_tree.get<double >("poissonratio.value",0.45);
  double YoungModulus=para_tree.get<double >("youngmodulus.value",1e5);
  dmetric=0.01; size_t length=12; size_t width=12; size_t height=12;

  printf("Poissonratio:: %lf  Youngmodulus:: %lf\n ",PoissonRatio,YoungModulus);
  size_t i,j,k;
  double material[2];
  double energy;
  material[0]=YoungModulus/(1+PoissonRatio)*0.5;
  material[1]=YoungModulus*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);
  for(i=0;i<3;++i)
    {
      for(j=i;j<3;++j)
	{
	  object* test_obj=new object(dmetric,length,width,height,100,1000,0,0,"linear");
	  configMaterial(test_obj,material);
	  myvector e1=myvector(0,0,0); myvector e2=myvector(0,0,0);
	  e1(i)=1; e2(j)=1;     
	  configForce(test_obj,e1,e2);
	  test_obj->harmonic_def_static();
	  energy=test_obj->calElasticEnergy();
	  printf("energy is :%lf \n ",energy);
	  stringstream ss;string i_str;ss<<(i+1); ss>>i_str;
	  stringstream mm;string j_str;mm<<(j+1); mm>>j_str;
	  string path_name=out_dir+"/test_obj_harmonic_def_"+i_str+"_"+j_str+".vtk";
	  myio.saveAsVTK(test_obj,path_name);
	  delete test_obj;
	}
    }
  return 0;
}

//input_dir & out_dir 的层数都是到对应的object文件夹
int harmonic_solver::solve(boost::property_tree::ptree &para_tree)
{
  io myio=io();
  string input_dir=para_tree.get<string>("input_dir.value");
  string out_dir=para_tree.get<string>("out_dir.value");
  string object_name=para_tree.get<string>("object_name.value");
  int level = para_tree.get<int>("level.value");
  size_t i,j,k;
  double energy;
  
  string input_fine_object=input_dir+"/"+object_name+".sub1.vtk";
  string input_fine_object_mat=input_dir+"/"+object_name+".sub1.mat";
  dmetric=myio.getDmetric(input_fine_object);
  for(i=0;i<3;++i)
    {
      for(j=i;j<3;++j)
	{
	  // harmonic solver中使用的是linear without stiffness tensor
	  object* test_obj=new object(input_fine_object,100,1000,0,0,"linear");
	  myio.getMatPara(test_obj,input_fine_object_mat,1);
	  myvector e1=myvector(0,0,0); myvector e2=myvector(0,0,0);

	  /*
	  if(i==j)
	    {
	      e1(i)=1; e2(j)=1;
	    }
	  else
	    {
	      e1(i)=-1; e2(j)=1;
	      }*/
	  e1(i)=1; e2(j)=1;
	  // e1(i)=-1; e2(j)=1;
	  configForce(test_obj,e1,e2);
	  test_obj->harmonic_def_static();
	  energy=test_obj->calElasticEnergy();
	  printf("energy is :%lf \n ",energy);
	  stringstream ss;string i_str;ss<<(i+1); ss>>i_str;
	  stringstream mm;string j_str;mm<<(j+1); mm>>j_str;
	  string path_name=out_dir+"/"+object_name+"_harmonic_def_"+i_str+"_"+j_str+".vtk";
	  myio.saveAsVTK(test_obj,path_name);
	  delete test_obj;
	}
    }
  return 0;
}

int harmonic_solver::configMaterial(object *myobject,const double (&material_para)[2])
{
  size_t num_hexahedrons=myobject->num_hexahedrons;
  size_t i,j,k,a,b,c;
  for(i=0;i<num_hexahedrons;++i)
    {
      for(a=0;a<2;++a)
	{
	  for(b=0;b<2;++b)
	    {
	      for(c=0;c<2;++c)
		{
		  myobject->myhexahedrons[i].material_para[a][b][c][0]=material_para[0];
		  myobject->myhexahedrons[i].material_para[a][b][c][1]=material_para[1];
		}
	    }
	}
    }
  return 0;
}

int harmonic_solver::configForce(object *myobject, myvector e1, myvector e2)
{
  size_t i,j,k;
  MatrixXd strain_harmonic(3,3);
  MatrixXd normal_matrix(3,1);
  MatrixXd force_matrix(3,1);
  for(i=0;i<3;++i)
    {
      for(j=0;j<3;++j)
	{
      	  strain_harmonic(i,j)=0.5*(e1(i)*e2(j)+e2(i)*e1(j));
	}
    }
  double area_face=pow(dmetric,2);
  for(i=0;i<myobject->num_faces;++i)
    {
      if(myobject->myfaces[i].num_hex==1)
	{
	  for(j=0;j<3;++j)
	    {
	      normal_matrix(j,0)=myobject->myfaces[i].normal_ori(j);
	    }
	  force_matrix=0.25*area_face*(strain_harmonic*normal_matrix);
	  for(j=0;j<4;++j)
	    {
	      myobject->myvertexs[myobject->myfaces[i].index_vertex[j]].force_external+=scale*myvector(force_matrix(0,0),force_matrix(1,0),force_matrix(2,0));
	    }
	}
    }
  return 0;
}
