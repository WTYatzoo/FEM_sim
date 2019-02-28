#include "io.h"
using namespace std;
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

int io::saveAsVTKwithForce(object* myobject,const string name)
{
  FILE *fp;
  size_t i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"tet\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",myobject->num_vertex);
  for(i=0;i<myobject->num_vertex;++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myobject->myvertexs[i].location.x,myobject->myvertexs[i].location.y,myobject->myvertexs[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",myobject->num_hexahedrons,myobject->num_hexahedrons*9);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"8");
      for(j=0;j<8;++j)
	{
	  fprintf(fp," %d",myobject->myhexahedrons[i].index_vertex[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"CELL_TYPES %d\n",myobject->num_hexahedrons);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"12\n");
    }
  fprintf(fp,"POINT_DATA %d\n",myobject->num_vertex);
  fprintf(fp,"VECTORS force double\n");
  //  fprintf(fp,"LOOKUP_TABLE default\n"); //writing this will be wrong!!!
  for(i=0;i<myobject->num_vertex;++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myobject->myvertexs[i].force_external.x,myobject->myvertexs[i].force_external.y,myobject->myvertexs[i].force_external.z);
    }
  fclose(fp);
  return 0;
}

int io::saveAsVTK(object* myobject,const string name)
{
  FILE *fp;
  size_t i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"tet\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",myobject->num_vertex);
  for(i=0;i<myobject->num_vertex;++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myobject->myvertexs[i].location.x,myobject->myvertexs[i].location.y,myobject->myvertexs[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",myobject->num_hexahedrons,myobject->num_hexahedrons*9);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"8");
      for(j=0;j<8;++j)
	{
	  fprintf(fp," %d",myobject->myhexahedrons[i].index_vertex[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"CELL_TYPES %d\n",myobject->num_hexahedrons);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"12\n");
    }
  
  fclose(fp);
  return 0;
}

int io::saveArrayAsVTK(std::vector<vertex> &myvertexs,std::vector<hexahedron> &myhexahedrons,const std::string name)
{
  FILE *fp;
  size_t i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"tet\n");
  fprintf(fp,"ASCII\n\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",myvertexs.size());
  for(i=0;i<myvertexs.size();++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myvertexs[i].location.x,myvertexs[i].location.y,myvertexs[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",myhexahedrons.size(),myhexahedrons.size()*9);
  for(i=0;i<myhexahedrons.size();++i)
    {
      fprintf(fp,"8");
      for(j=0;j<8;++j)
	{
	  fprintf(fp," %d",myhexahedrons[i].index_vertex[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"CELL_TYPES %d\n",myhexahedrons.size());
  for(i=0;i<myhexahedrons.size();++i)
    {
      fprintf(fp,"12\n");
    }
  fclose(fp);
  return 0;
}

int io::saveAsVTKwithPara(object* myobject,const string name)
{
  FILE *fp;
  size_t i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"tet\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",myobject->num_vertex);
  for(i=0;i<myobject->num_vertex;++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myobject->myvertexs[i].location.x,myobject->myvertexs[i].location.y,myobject->myvertexs[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",myobject->num_hexahedrons,myobject->num_hexahedrons*9);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"8");
      for(j=0;j<8;++j)
	{
	  fprintf(fp," %d",myobject->myhexahedrons[i].index_vertex[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"CELL_TYPES %d\n",myobject->num_hexahedrons);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"12\n");
    }

  fprintf(fp,"CELL_DATA %d\n",myobject->num_hexahedrons);
  fprintf(fp,"SCALARS para double\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"%lf\n",myobject->myhexahedrons[i].material_para[0][0][0][0]);
    }
  fclose(fp);
  return 0;
}

int io::saveAsVTKwithScalarPara(object* myobject,int which,const string name)
{
  FILE *fp;
  size_t i,j;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"tet\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d double\n",myobject->num_vertex);
  for(i=0;i<myobject->num_vertex;++i)
    {
      fprintf(fp,"%lf %lf %lf\n",myobject->myvertexs[i].location.x,myobject->myvertexs[i].location.y,myobject->myvertexs[i].location.z);
    }
  fprintf(fp,"CELLS %d %d\n",myobject->num_hexahedrons,myobject->num_hexahedrons*9);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"8");
      for(j=0;j<8;++j)
	{
	  fprintf(fp," %d",myobject->myhexahedrons[i].index_vertex[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]]);
	}
      fprintf(fp,"\n");

    }
  fprintf(fp,"CELL_TYPES %d\n",myobject->num_hexahedrons);
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"12\n");
    }

  fprintf(fp,"CELL_DATA %d\n",myobject->num_hexahedrons);
  fprintf(fp,"SCALARS para double\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(i=0;i<myobject->num_hexahedrons;++i)
    {
      fprintf(fp,"%lf\n",myobject->myhexahedrons[i].ScalarPara(which,0));
    }
  fclose(fp);
  return 0;
}

int io::getMatPara(object* myobject,const std::string name,int kind)
{
  FILE *fp;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"r");
  size_t i;
  size_t a,b,c,d;
  size_t num_hexahedrons=myobject->num_hexahedrons;
  printf("num_hexahedrons: %d\n",myobject->num_hexahedrons);
  if(kind==1)
    {
      for(i=0;i<num_hexahedrons;++i)
	{
	  for(a=0;a<2;++a)
	    {
	      for(b=0;b<2;++b)
		{
		  for(c=0;c<2;++c)
		    {
		      for(d=0;d<2;++d)
			{
			  fscanf(fp,"%lf ",&myobject->myhexahedrons[i].material_para[a][b][c][d]);
			}
		    }
		}
	    }
	  fscanf(fp,"\n");
	}
    }
  else if(kind==0)
    {
      for(i=0;i<num_hexahedrons;++i)
	{
	  for(a=0;a<6;++a)
	    {
	      for(b=0;b<6;++b)
		{
		  fscanf(fp,"%lf ",&myobject->myhexahedrons[i].stiffness_tensor(a,b));
		}
	    }
	  fscanf(fp,"\n");
	}
    }
  fclose(fp);
  return 0;
}
int io::saveMatPara(object* myobject,const string name,int kind)
{
  FILE *fp;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  size_t i;
  size_t a,b,c,d;
  size_t num_hexahedrons=myobject->num_hexahedrons;
  if(kind==1)
    {
      for(i=0;i<num_hexahedrons;++i)
	{
	  for(a=0;a<2;++a)
	    {
	      for(b=0;b<2;++b)
		{
		  for(c=0;c<2;++c)
		    {
		      for(d=0;d<2;++d)
			{
			  fprintf(fp,"%lf ",myobject->myhexahedrons[i].material_para[a][b][c][d]);
			}
		    }
		}
	    }
	  fprintf(fp,"\n");
	}
    }
  else if(kind==0)
    {
      for(i=0;i<num_hexahedrons;++i)
	{
	  for(a=0;a<6;++a)
	    {
	      for(b=0;b<6;++b)
		{
		  fprintf(fp,"%lf ",myobject->myhexahedrons[i].stiffness_tensor(a,b));
		}
	    }
	  fprintf(fp,"\n");
	}
    }
  fclose(fp);
  return 0;
}

int io::getConstraintFromCsv(object* myobject,const std::string name)
{
  FILE *fp;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"r");
  char filter[400];
  fscanf(fp,"%s",filter);
  double x[3];
  int fixed_now;
  myobject->num_fixed=0;
  
  while(!feof(fp))
    {
      fscanf(fp,"%d,%lf,%lf,%lf\n",&fixed_now,&x[0],&x[1],&x[2]);
      //      printf("fixed: %d\n",fixed_now);
      myobject->myvertexs[fixed_now].isFixed=1;
      myobject->num_fixed++;
    }
  fclose(fp);
  return 0;
}

double io::getDmetric(const std::string name)
{
  printf("%s \n",name.c_str());
  /*
  FILE* fp;
  size_t i,j;
  fp=fopen(name.c_str(),"r");
  char filter[10];
  do
   {
     fscanf(fp,"%s",filter);
   }while(filter[0]!='P');
  
  int num_vertex_file;
  fscanf(fp,"%d%s",&num_vertex_file,filter);
  myvector points[2];
  for(i=0;i<2;++i)
    {
      for(j=0;j<3;++j)
	{
	  fscanf(fp,"%lf",&(points[i](j)));
	}
    }
  fclose(fp);
  return (points[0]-points[1]).len();*/
  object* test_obj=new object();
  getVertexAndHex(test_obj->myvertexs,test_obj->myhexahedrons,name);
  myvector loc[2];
  loc[0]=test_obj->myvertexs[test_obj->myhexahedrons[0].index_vertex[0][0][0]].location_original;
  loc[1]=test_obj->myvertexs[test_obj->myhexahedrons[0].index_vertex[1][0][0]].location_original;
  double dmetric=(loc[0]-loc[1]).len();
  delete test_obj;
  return dmetric;
}

int io::getVertexAndHex(std::vector<vertex> &myvertexs,std::vector<hexahedron> &myhexahedrons,const std::string name)
{
  printf("%s \n",name.c_str());
  FILE* fp;
  fp=fopen(name.c_str(),"r");
  char filter[10];
  do
   {
     fscanf(fp,"%s",filter);
   }while(filter[0]!='P');

  int num_vertex_file;
  int index_vertex_file[2][2][2]; size_t index_vertex[2][2][2];
  int num_hexahedrons_file,num_hexahedrons9_file;
  int filter8;
  // 从文件中读取整数用int接受%d ,即使使用时是用size_t ,如果直接用size_t接受%u使用是会出错，此处可能是本身的bug
  double x,y,z;

  fscanf(fp,"%d%s",&num_vertex_file,filter);
  size_t i,j;
  for(i=0;i<num_vertex_file;++i)
    {
      fscanf(fp,"%lf%lf%lf",&x,&y,&z);
      myvertexs.push_back(vertex(myvector(x,y,z)));
    }
  fscanf(fp,"%s%d%d",filter,&num_hexahedrons_file,&num_hexahedrons9_file);
 
  for(i=0;i<num_hexahedrons_file;++i)
    {
      fscanf(fp,"%d",&filter8);
      for(j=0;j<8;++j)
	{
	  fscanf(fp,"%d",&index_vertex_file[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]]);
	  index_vertex[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]]=index_vertex_file[vertex_index[j][0]][vertex_index[j][1]][vertex_index[j][2]];
	}
      myhexahedrons.push_back(hexahedron(index_vertex));
    }
  fclose(fp);
  return 0;
}


int io::getMatParaFromTXT(object* myobject,const std::string name)
{
  FILE* file=fopen(name.c_str(),"r");
  int which_hex; double y,p;
  size_t a,b,c,d;
  while(!feof(file))
    {
      fscanf(file,"%d %lf %lf",&which_hex,&y,&p);
      double lame1=y/(1+p)*0.5;
      double lame2=y*p/(1+p)/(1-2*p);
      for(a=0;a<2;++a)
	{
	  for(b=0;b<2;++b)
	    {
	      for(c=0;c<2;++c)
		{
		  myobject->myhexahedrons[which_hex].material_para[a][b][c][0]=lame1;
		  myobject->myhexahedrons[which_hex].material_para[a][b][c][1]=lame2;
		}
	    }
	}
    }
  fclose(file);
  return 0;
}

int io::saveMatParaAsTXT(object* myobject,const std::string name)
{
  FILE* file=fopen(name.c_str(),"w");
  double y,p,lame1,lame2;
  size_t i;
  for(i=0;i<myobject->num_hexahedrons;i++)
    {
      lame1=myobject->myhexahedrons[i].material_para[0][0][0][0];
      lame2=myobject->myhexahedrons[i].material_para[0][0][0][1];
      p=lame2/(lame1+lame2)*0.5;
      y=2*lame1*(p+1);
      fprintf(file,"%d %lf %lf\n",i,y,p);
    }
  fclose(file);
  return 0;
}

int io::save_ERROR_TXT(double error_value,const std::string name)
{
  FILE* file=fopen(name.c_str(),"a");
  fprintf(file,"%lf\n",error_value);
  fclose(file);
  return 0;
}

int io::saveTimeNormPair(object* myobject,const std::string name)
{
  FILE *fp;
  printf("%s \n",name.c_str());
  fp=fopen(name.c_str(),"w");
  // fprintf(fp,"time,norm_Jacobian_cal\n");
  vector<pair<double,double > >::iterator it;
  for(it=myobject->time_norm_pair.begin();it!=myobject->time_norm_pair.end();it++)
    {
      fprintf(fp,"%lf %lf\n",it->first,it->second);
    }
  fclose(fp);
  return 0;
}

