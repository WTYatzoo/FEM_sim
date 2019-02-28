#include "object_creator.h"
#include "io.h"
using namespace std;
using namespace Eigen;

class R
{
public:
  size_t index_hex;
  double x_center;
  double y,p;
  ~R(){}
  R(){}
  R(size_t index_hex,double x_center)
    {
      this->index_hex=index_hex;
      this->x_center=x_center;
    }
};

bool cmp1(const R&r1,const R&r2)
{
    return r1.x_center<r2.x_center;
}
bool cmp2(const R&r1,const R&r2)
{
    return r1.index_hex<r2.index_hex;
}


object_creator::object_creator(boost::property_tree::ptree &para_tree)
{
  out_dir=para_tree.get<string>("out_dir.value");
  object_name=para_tree.get<string>("object_name.value");
  dmetric_fine=para_tree.get<double>("dmetric_fine.value");

  level=para_tree.get<int>("level.value");
  //此时lc_fine wc_fine hc_fine 不表示节点数
  lc_fine=para_tree.get<int>("l_size_fine.value");
  wc_fine=para_tree.get<int>("w_size_fine.value");
  hc_fine=para_tree.get<int>("h_size_fine.value");
  
  PoissonRatio=para_tree.get<double >("poissonratio.value",0.45);
  YoungModulus1=para_tree.get<double >("youngmodulus1.value",1e4); //hard 
  YoungModulus2=para_tree.get<double >("youngmodulus2.value",1e3); //soft

  //不涉及 constitutive_model
  //create_fine_object();
  create_object();
  // create_object_bridge();
  //create_object_bridge_sin();
}

object_creator::object_creator()
{
  ;
}
object_creator::~object_creator()
{
   if(index_for_vertex!=NULL)
    {
      size_t i,j;
      for(i=0;i<lc_fine;++i)
	{
	  for(j=0;j<wc_fine;++j)
	    {
	      delete[] index_for_vertex[i][j];
	    }
	  delete[] index_for_vertex[i];
	}
      delete[] index_for_vertex;
    }
}

//用来做中间软两边硬的两端固定住的棒子 材料分布为y=x^2
int object_creator::create_object_bridge()
{
  io myio=io();
  object* fine_obj=new object();
  object* coarsen_obj=new object();

  dmetric_coarsen=dmetric_fine*level;
  length_fine=dmetric_fine*lc_fine; width_fine=dmetric_fine*wc_fine; height_fine=dmetric_fine*hc_fine;
  lc_fine++; wc_fine++; hc_fine++;
  
  this->index_for_vertex=NULL;
  double length_now,width_now,height_now;
  size_t i,j,k;
  size_t x,y,z;
  size_t a,b,c,d;

  // 分配内存空间 保存顶点索引
  index_for_vertex=(int***)new int**[lc_fine];
  for(i=0;i<lc_fine;++i)
    {
      index_for_vertex[i]=(int**)new int*[wc_fine];
      for(j=0;j<wc_fine;++j)
	{
	  index_for_vertex[i][j]=new int[hc_fine];
	}
    }
  
  for(i=0;i<lc_fine;++i)
    {
      for(j=0;j<wc_fine;++j)
	{
	  for(k=0;k<hc_fine;++k)
	    {
	      index_for_vertex[i][j][k]=-1;
	    }
	}
    }
  
  size_t index_now_vertex=0;
  //set index for all vertexs and EPS is to make double`error correction and set it as local constant
  double EPS=1e-6;
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_coarsen,i+=level)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_coarsen,j+=level)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_coarsen,k+=level)
	    {
	      index_for_vertex[i][j][k]=index_now_vertex;
	      coarsen_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	      index_now_vertex++;
	    }
	}
    }

  coarsen_obj->num_vertex=coarsen_obj->myvertexs.size();
  for(i=0;i<coarsen_obj->num_vertex;++i)
    {
      fine_obj->myvertexs.push_back(coarsen_obj->myvertexs[i]);
    }
  
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_fine,++i)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_fine,++j)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_fine,++k)
	    {
	      if(index_for_vertex[i][j][k]==-1)
		{
		  index_for_vertex[i][j][k]=index_now_vertex;
		  fine_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	          index_now_vertex++;
		}	      	      
	    }
	}
    }   
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  
  size_t index_vertex_now[2][2][2];
  //traverse all the grid

  double material_para[2][2][2][2];
  for(i=0;i<lc_fine-level;i+=level)
    {   
      for(j=0;j<wc_fine-level;j+=level)
	{ 
	  for(k=0;k<hc_fine-level;k+=level)
	    {
	      for(x=0;x<level;++x)
		{
		  for(y=0;y<level;++y)
		    {
		      for(z=0;z<level;++z)
			{
			  for(a=0;a<2;++a)
			    {
			      for(b=0;b<2;++b)
				{
				  for(c=0;c<2;++c)
				    {
				      index_vertex_now[a][b][c]=index_for_vertex[i+x+a][j+y+b][k+z+c];
				    }
				}
			    }
			  fine_obj->myhexahedrons.push_back(hexahedron(index_vertex_now));
			}
		    }
		}
	   
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  index_vertex_now[a][b][c]=index_for_vertex[i+a*level][j+b*level][k+c*level];
			}
		    }
		}
	      coarsen_obj->myhexahedrons.push_back(hexahedron(index_vertex_now));
	    }    
	}
    }
  
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  coarsen_obj->num_hexahedrons=coarsen_obj->myhexahedrons.size();

  { 
    vector< R > my_R;
    while(!my_R.empty())
      {
	my_R.pop_back();
      }
    for(i=0;i<fine_obj->num_hexahedrons;i++)
      {
	R r_now=R(0,0);
	double x_center_now=0.0; 
	for(a=0;a<2;a++)
	  {
	    for(b=0;b<2;b++)
	      {
		for(c=0;c<2;c++)
		  {
		    x_center_now+=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[a][b][c]].location(0);
		  }
	      }
	  }
	x_center_now*=0.125;
	r_now.x_center=x_center_now; r_now.index_hex=i;
	my_R.push_back(r_now);
      }
    int len=my_R.size();
    sort(my_R.begin(),my_R.end(),cmp1);

    double x_center_min=my_R[0].x_center;
    double x_center_max=my_R[my_R.size()-1].x_center;
    double x_center_middle=0.5*(x_center_max+x_center_min);
    double dif_min=fabs(x_center_min-x_center_middle);
    double x_center_middle_actual;
    for(i=0;i<my_R.size();i++)
      {
	if(fabs(my_R[i].x_center-x_center_middle)<=dif_min)
	  {
	    x_center_middle_actual=my_R[i].x_center; 
	    dif_min=fabs(my_R[i].x_center-x_center_middle);
	  }
	else
	  {
	    ;
	  }
      }
    double EPS=1e-5;
    for(i=0;i<my_R.size();i++)
      {
	if(fabs(my_R[i].x_center-x_center_middle_actual)<=EPS)
	  {
	    my_R[i].y=YoungModulus2; my_R[i].p=PoissonRatio;
	  }
	else if(fabs(my_R[i].x_center-x_center_min)<=EPS)
	  {
	    my_R[i].y=YoungModulus1; my_R[i].p=PoissonRatio;
	  }
	else if(fabs(my_R[i].x_center-x_center_max)<=EPS)
	  {
	    my_R[i].y=YoungModulus1; my_R[i].p=PoissonRatio;
	  }
	else
	  {
	    my_R[i].y=-1; my_R[i].p=PoissonRatio;
	  }
      }

    double kk=(YoungModulus1-YoungModulus2)/pow((x_center_max-x_center_middle),2);
    for(i=0;i<my_R.size();i++)
      {
	if(my_R[i].y>0.0)
	  {
	    ;
	  }
	else
	  {
	    my_R[i].y=kk*pow(my_R[i].x_center-x_center_middle,2)+YoungModulus2;
	  }
      }
    sort(my_R.begin(),my_R.end(),cmp2);

    for(i=0;i<len;i++)
      {
	for(a=0;a<2;a++)
	  {
	    for(b=0;b<2;b++)
	      {
		for(c=0;c<2;c++)
		  {
		    fine_obj->myhexahedrons[i].material_para[a][b][c][0]=my_R[i].y/(1+PoissonRatio)*0.5;
		    fine_obj->myhexahedrons[i].material_para[a][b][c][1]=my_R[i].y*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);
		  }
	      }
	  }
      }
  }

  size_t ct_fine_hex=0,ct_coarsen_hex=0;
  for(i=0;i<lc_fine-level;i+=level)
    {   
      for(j=0;j<wc_fine-level;j+=level)
	{ 
	  for(k=0;k<hc_fine-level;k+=level)
	    {
	      for(x=0;x<level;++x)
		{
		  for(y=0;y<level;++y)
		    {
		      for(z=0;z<level;++z)
			{
			  for(a=0;a<2;++a)
			    {
			      for(b=0;b<2;++b)
				{
				  for(c=0;c<2;++c)
				    {
				      index_vertex_now[a][b][c]=index_for_vertex[i+x+a][j+y+b][k+z+c];
				    }
				}
			    }
			  material_para[x][y][z][0]=fine_obj->myhexahedrons[ct_fine_hex].material_para[0][0][0][0];
			  material_para[x][y][z][1]=fine_obj->myhexahedrons[ct_fine_hex].material_para[0][0][0][1];
			  ct_fine_hex++;
			}
		    }
		}
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  index_vertex_now[a][b][c]=index_for_vertex[i+a*level][j+b*level][k+c*level];
			}
		    }
		}
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  for(d=0;d<2;++d)
			    {
			      coarsen_obj->myhexahedrons[ct_coarsen_hex].material_para[a][b][c][d]=material_para[a][b][c][d];
			    }			  
			}
		    }
		}
	      ct_coarsen_hex++;
	    }    
	}
    }
  i=1;
  stringstream ss_fine; string sub_fine; ss_fine<<i; ss_fine>>sub_fine;
  string path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".vtk";
  myio.saveAsVTKwithPara(fine_obj,path_name_fine);
  string mat_path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".mat";
  myio.saveMatPara(fine_obj,mat_path_name_fine,1);// 1表示save material_para array ; 0表示save stiffness_tensor;

  i=0;
  stringstream ss_coarsen; string sub_coarsen; ss_coarsen<<i; ss_coarsen>>sub_coarsen;
  string path_name_coarsen=out_dir+"/"+object_name+".sub"+sub_coarsen+".vtk";
  myio.saveAsVTKwithPara(coarsen_obj,path_name_coarsen);
  string mat_path_name_coarsen=out_dir+"/"+object_name+".sub"+sub_coarsen+".mat"; //1.简单平均得到的粗网格材质 或者　2.直接拿过来
  myio.saveMatPara(coarsen_obj,mat_path_name_coarsen,1);// 1表示save material_para array ; 0表示save stiffness_tensor;
  
  delete fine_obj;
  delete coarsen_obj;
  
  return 0;
}

//用来做中间软两边硬的两端固定住的棒子 材料分布符合y=sin(x)
int object_creator::create_object_bridge_sin()
{
  io myio=io();
  object* fine_obj=new object();
  object* coarsen_obj=new object();

  dmetric_coarsen=dmetric_fine*level;
  length_fine=dmetric_fine*lc_fine; width_fine=dmetric_fine*wc_fine; height_fine=dmetric_fine*hc_fine;
  lc_fine++; wc_fine++; hc_fine++;
  
  this->index_for_vertex=NULL;
  double length_now,width_now,height_now;
  size_t i,j,k;
  size_t x,y,z;
  size_t a,b,c,d;

  // 分配内存空间 保存顶点索引
  index_for_vertex=(int***)new int**[lc_fine];
  for(i=0;i<lc_fine;++i)
    {
      index_for_vertex[i]=(int**)new int*[wc_fine];
      for(j=0;j<wc_fine;++j)
	{
	  index_for_vertex[i][j]=new int[hc_fine];
	}
    }
  
  for(i=0;i<lc_fine;++i)
    {
      for(j=0;j<wc_fine;++j)
	{
	  for(k=0;k<hc_fine;++k)
	    {
	      index_for_vertex[i][j][k]=-1;
	    }
	}
    }
  
  size_t index_now_vertex=0;
  //set index for all vertexs and EPS is to make double`error correction and set it as local constant
  double EPS=1e-6;
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_coarsen,i+=level)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_coarsen,j+=level)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_coarsen,k+=level)
	    {
	      index_for_vertex[i][j][k]=index_now_vertex;
	      coarsen_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	      index_now_vertex++;
	    }
	}
    }

  coarsen_obj->num_vertex=coarsen_obj->myvertexs.size();
  for(i=0;i<coarsen_obj->num_vertex;++i)
    {
      fine_obj->myvertexs.push_back(coarsen_obj->myvertexs[i]);
    }
  
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_fine,++i)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_fine,++j)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_fine,++k)
	    {
	      if(index_for_vertex[i][j][k]==-1)
		{
		  index_for_vertex[i][j][k]=index_now_vertex;
		  fine_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	          index_now_vertex++;
		}	      	      
	    }
	}
    }   
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  
  size_t index_vertex_now[2][2][2];
  //traverse all the grid

  double material_para[2][2][2][2];
  for(i=0;i<lc_fine-level;i+=level)
    {   
      for(j=0;j<wc_fine-level;j+=level)
	{ 
	  for(k=0;k<hc_fine-level;k+=level)
	    {
	      for(x=0;x<level;++x)
		{
		  for(y=0;y<level;++y)
		    {
		      for(z=0;z<level;++z)
			{
			  for(a=0;a<2;++a)
			    {
			      for(b=0;b<2;++b)
				{
				  for(c=0;c<2;++c)
				    {
				      index_vertex_now[a][b][c]=index_for_vertex[i+x+a][j+y+b][k+z+c];
				    }
				}
			    }
			  fine_obj->myhexahedrons.push_back(hexahedron(index_vertex_now));
			}
		    }
		}
	   
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  index_vertex_now[a][b][c]=index_for_vertex[i+a*level][j+b*level][k+c*level];
			}
		    }
		}
	      coarsen_obj->myhexahedrons.push_back(hexahedron(index_vertex_now));
	    }    
	}
    }
  
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  coarsen_obj->num_hexahedrons=coarsen_obj->myhexahedrons.size();

  { 
    vector< R > my_R;
    while(!my_R.empty())
      {
	my_R.pop_back();
      }
    for(i=0;i<fine_obj->num_hexahedrons;i++)
      {
	R r_now=R(0,0);
	double x_center_now=0.0; 
	for(a=0;a<2;a++)
	  {
	    for(b=0;b<2;b++)
	      {
		for(c=0;c<2;c++)
		  {
		    x_center_now+=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[a][b][c]].location(0);
		  }
	      }
	  }
	x_center_now*=0.125;
	r_now.x_center=x_center_now; r_now.index_hex=i;
	my_R.push_back(r_now);
      }
    int len=my_R.size();
    sort(my_R.begin(),my_R.end(),cmp1);

    double x_center_min=my_R[0].x_center;
    double x_center_max=my_R[my_R.size()-1].x_center;
    double EPS=1e-5;

    for(i=0;i<my_R.size();i++)
      {
	my_R[i].y=(YoungModulus1-YoungModulus2)*0.5*sin(50/(x_center_max-x_center_min)*my_R[i].x_center-50*x_center_min/(x_center_max-x_center_min))+(YoungModulus2+YoungModulus1)*0.5; // 2*pai*8约等于50
	my_R[i].p=PoissonRatio;
      }
    sort(my_R.begin(),my_R.end(),cmp2);

    for(i=0;i<len;i++)
      {
	for(a=0;a<2;a++)
	  {
	    for(b=0;b<2;b++)
	      {
		for(c=0;c<2;c++)
		  {
		    fine_obj->myhexahedrons[i].material_para[a][b][c][0]=my_R[i].y/(1+PoissonRatio)*0.5;
		    fine_obj->myhexahedrons[i].material_para[a][b][c][1]=my_R[i].y*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);
		  }
	      }
	  }
      }
  }

  size_t ct_fine_hex=0,ct_coarsen_hex=0;
  for(i=0;i<lc_fine-level;i+=level)
    {   
      for(j=0;j<wc_fine-level;j+=level)
	{ 
	  for(k=0;k<hc_fine-level;k+=level)
	    {
	      for(x=0;x<level;++x)
		{
		  for(y=0;y<level;++y)
		    {
		      for(z=0;z<level;++z)
			{
			  for(a=0;a<2;++a)
			    {
			      for(b=0;b<2;++b)
				{
				  for(c=0;c<2;++c)
				    {
				      index_vertex_now[a][b][c]=index_for_vertex[i+x+a][j+y+b][k+z+c];
				    }
				}
			    }
			  material_para[x][y][z][0]=fine_obj->myhexahedrons[ct_fine_hex].material_para[0][0][0][0];
			  material_para[x][y][z][1]=fine_obj->myhexahedrons[ct_fine_hex].material_para[0][0][0][1];
			  ct_fine_hex++;
			}
		    }
		}
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  index_vertex_now[a][b][c]=index_for_vertex[i+a*level][j+b*level][k+c*level];
			}
		    }
		}
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  for(d=0;d<2;++d)
			    {
			      coarsen_obj->myhexahedrons[ct_coarsen_hex].material_para[a][b][c][d]=material_para[a][b][c][d];
			    }			  
			}
		    }
		}
	      ct_coarsen_hex++;
	    }    
	}
    }
  i=1;
  stringstream ss_fine; string sub_fine; ss_fine<<i; ss_fine>>sub_fine;
  string path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".vtk";
  myio.saveAsVTKwithPara(fine_obj,path_name_fine);
  string mat_path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".mat";
  myio.saveMatPara(fine_obj,mat_path_name_fine,1);// 1表示save material_para array ; 0表示save stiffness_tensor;

  i=0;
  stringstream ss_coarsen; string sub_coarsen; ss_coarsen<<i; ss_coarsen>>sub_coarsen;
  string path_name_coarsen=out_dir+"/"+object_name+".sub"+sub_coarsen+".vtk";
  myio.saveAsVTKwithPara(coarsen_obj,path_name_coarsen);
  string mat_path_name_coarsen=out_dir+"/"+object_name+".sub"+sub_coarsen+".mat"; //1.简单平均得到的粗网格材质 或者　2.直接拿过来
  myio.saveMatPara(coarsen_obj,mat_path_name_coarsen,1);// 1表示save material_para array ; 0表示save stiffness_tensor;
  
  delete fine_obj;
  delete coarsen_obj;
  
  return 0;
}

int object_creator::create_object()
{
  io myio=io();
  object* fine_obj=new object();
  object* coarsen_obj=new object();

  dmetric_coarsen=dmetric_fine*level;
  length_fine=dmetric_fine*lc_fine; width_fine=dmetric_fine*wc_fine; height_fine=dmetric_fine*hc_fine;
  lc_fine++; wc_fine++; hc_fine++;
  
  this->index_for_vertex=NULL;
  double length_now,width_now,height_now;
  size_t i,j,k;
  size_t x,y,z;
  size_t a,b,c,d;

  // 分配内存空间 保存顶点索引
  index_for_vertex=(int***)new int**[lc_fine];
  for(i=0;i<lc_fine;++i)
    {
      index_for_vertex[i]=(int**)new int*[wc_fine];
      for(j=0;j<wc_fine;++j)
	{
	  index_for_vertex[i][j]=new int[hc_fine];
	}
    }
  
  for(i=0;i<lc_fine;++i)
    {
      for(j=0;j<wc_fine;++j)
	{
	  for(k=0;k<hc_fine;++k)
	    {
	      index_for_vertex[i][j][k]=-1;
	    }
	}
    }
  
  size_t index_now_vertex=0;
  //set index for all vertexs and EPS is to make double`error correction and set it as local constant
  double EPS=1e-6;
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_coarsen,i+=level)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_coarsen,j+=level)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_coarsen,k+=level)
	    {
	      index_for_vertex[i][j][k]=index_now_vertex;
	      coarsen_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	      index_now_vertex++;
	    }
	}
    }

  coarsen_obj->num_vertex=coarsen_obj->myvertexs.size();
  for(i=0;i<coarsen_obj->num_vertex;++i)
    {
      fine_obj->myvertexs.push_back(coarsen_obj->myvertexs[i]);
    }
  
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_fine,++i)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_fine,++j)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_fine,++k)
	    {
	      if(index_for_vertex[i][j][k]==-1)
		{
		  index_for_vertex[i][j][k]=index_now_vertex;
		  fine_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	          index_now_vertex++;
		}	      	      
	    }
	}
    }   
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  
  size_t index_vertex_now[2][2][2];
  //traverse all the grid

  double YoungModulus_blend;

  double material_para[2][2][2][2];
  for(i=0;i<lc_fine-level;i+=level)
    {   
      for(j=0;j<wc_fine-level;j+=level)
	{ 
	  for(k=0;k<hc_fine-level;k+=level)
	    {
	      YoungModulus_blend=0;
	      for(x=0;x<level;++x)
		{
		  if(x==0)
		    {
		      YoungModulus_blend+=YoungModulus1;
		    }
		  else if(x==1&&i!=lc_fine-3)
		    {
		      YoungModulus_blend+=YoungModulus2;
		    }
		  else if(x==1&&i==lc_fine-3)
		    {
		      YoungModulus_blend+=YoungModulus1;
		    }
		  for(y=0;y<level;++y)
		    {
		      for(z=0;z<level;++z)
			{
			  for(a=0;a<2;++a)
			    {
			      for(b=0;b<2;++b)
				{
				  for(c=0;c<2;++c)
				    {
				      index_vertex_now[a][b][c]=index_for_vertex[i+x+a][j+y+b][k+z+c];
				    }
				}
			    }
			  fine_obj->myhexahedrons.push_back(hexahedron(index_vertex_now));
			  if(x%2==0)
			    {
			      bind_mat(fine_obj->myhexahedrons[fine_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus1);
			      material_para[x][y][z][0]=YoungModulus1/(1+PoissonRatio)*0.5;
			      material_para[x][y][z][1]=YoungModulus1*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);
			    }
			  else if(x%2==1&&i!=lc_fine-3)
			    {
			      bind_mat(fine_obj->myhexahedrons[fine_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus2);
			      material_para[x][y][z][0]=YoungModulus2/(1+PoissonRatio)*0.5;
			      material_para[x][y][z][1]=YoungModulus2*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);
			    }
			  else if(x%2==1&&i==lc_fine-3)
			    {
			      //就让它是ABABAB
			      bind_mat(fine_obj->myhexahedrons[fine_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus2);
			      material_para[x][y][z][0]=YoungModulus2/(1+PoissonRatio)*0.5;
			      material_para[x][y][z][1]=YoungModulus2*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);
			      // trick for 左右两边都是硬的　ABABAA
			      /*
			      bind_mat(fine_obj->myhexahedrons[fine_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus1);
			      material_para[x][y][z][0]=YoungModulus1/(1+PoissonRatio)*0.5;
			      material_para[x][y][z][1]=YoungModulus1*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);*/
			    }
			}
		    }
		}
	      YoungModulus_blend*=0.5;
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  index_vertex_now[a][b][c]=index_for_vertex[i+a*level][j+b*level][k+c*level];
			}
		    }
		}
	      coarsen_obj->myhexahedrons.push_back(hexahedron(index_vertex_now));

	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  for(d=0;d<2;++d)
			    {
			      coarsen_obj->myhexahedrons[coarsen_obj->myhexahedrons.size()-1].material_para[a][b][c][d]=material_para[a][b][c][d];
			    }			  
			}
		    }
		}
	      //    bind_mat(coarsen_obj->myhexahedrons[coarsen_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus_blend);
	    }    
	}
    }
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  coarsen_obj->num_hexahedrons=coarsen_obj->myhexahedrons.size();
  i=1;
  stringstream ss_fine; string sub_fine; ss_fine<<i; ss_fine>>sub_fine;
  string path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".vtk";
  myio.saveAsVTKwithPara(fine_obj,path_name_fine);
  string mat_path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".mat";
  myio.saveMatPara(fine_obj,mat_path_name_fine,1);// 1表示save material_para array ; 0表示save stiffness_tensor;

  i=0;
  stringstream ss_coarsen; string sub_coarsen; ss_coarsen<<i; ss_coarsen>>sub_coarsen;
  string path_name_coarsen=out_dir+"/"+object_name+".sub"+sub_coarsen+".vtk";
  myio.saveAsVTKwithPara(coarsen_obj,path_name_coarsen);
  string mat_path_name_coarsen=out_dir+"/"+object_name+".sub"+sub_coarsen+".mat"; //1.简单平均得到的粗网格材质 或者　2.直接拿过来
  myio.saveMatPara(coarsen_obj,mat_path_name_coarsen,1);// 1表示save material_para array ; 0表示save stiffness_tensor;
  
  delete fine_obj;
  delete coarsen_obj;
  
  return 0;
}

//暂时不会用了
int object_creator::create_fine_object()
{
  io myio=io();
  object* fine_obj=new object();
  dmetric_coarsen=dmetric_fine*level;
  length_fine=dmetric_fine*lc_fine; width_fine=dmetric_fine*wc_fine; height_fine=dmetric_fine*hc_fine;
  lc_fine++; wc_fine++; hc_fine++;
  
  this->index_for_vertex=NULL;
  double length_now,width_now,height_now;
  size_t i,j,k;
  size_t a,b,c;

  // 分配内存空间 保存顶点索引
  index_for_vertex=(int***)new int**[lc_fine];
  for(i=0;i<lc_fine;++i)
    {
      index_for_vertex[i]=(int**)new int*[wc_fine];
      for(j=0;j<wc_fine;++j)
	{
	  index_for_vertex[i][j]=new int[hc_fine];
	}
    }


  for(i=0;i<lc_fine;++i)
    {
      for(j=0;j<wc_fine;++j)
	{
	  for(k=0;k<hc_fine;++k)
	    {
	      index_for_vertex[i][j][k]=-1;
	    }
	}
    }
  
  size_t index_now_vertex=0;
  //set index for all vertexs and EPS is to make double`error correction and set it as local constant
  double EPS=1e-6;
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_coarsen,i+=level)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_coarsen,j+=level)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_coarsen,k+=level)
	    {
	      index_for_vertex[i][j][k]=index_now_vertex;
	      fine_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	      index_now_vertex++;
	    }
	}
    }
  
  for(i=0,length_now=-length_fine*0.5;length_now<length_fine*0.5+EPS;length_now+=dmetric_fine,++i)
    {
      for(j=0,width_now=-width_fine*0.5;width_now<width_fine*0.5+EPS;width_now+=dmetric_fine,++j)
	{
	  for(k=0,height_now=-height_fine*0.5;height_now<height_fine*0.5+EPS;height_now+=dmetric_fine,++k)
	    {
	      if(index_for_vertex[i][j][k]==-1)
		{
		  index_for_vertex[i][j][k]=index_now_vertex;
		  fine_obj->myvertexs.push_back(vertex(myvector(length_now,width_now,height_now)));
	          index_now_vertex++;
		}	      	      
	    }
	}
    }   
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  
  size_t index_vertex_now[2][2][2];
  //traverse all the grid
  
  for(i=0;i<lc_fine-1;++i)
    {   
      for(j=0;j<wc_fine-1;++j)
	{
	  for(k=0;k<hc_fine-1;++k)
	    {
	      for(a=0;a<2;++a)
		{
		  for(b=0;b<2;++b)
		    {
		      for(c=0;c<2;++c)
			{
			  index_vertex_now[a][b][c]=index_for_vertex[i+a][j+b][k+c];
			}
		    }
		}
	      fine_obj->myhexahedrons.push_back(hexahedron(index_vertex_now));
	      if(i%2==0)
		{
		  bind_mat(fine_obj->myhexahedrons[fine_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus1);
		}
	      else if(i%2==1&&((i+1)!=lc_fine-1))
		{
		  bind_mat(fine_obj->myhexahedrons[fine_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus2);
		}
	      else if(i%2==1&&((i+1)==lc_fine-1))
		{
		  bind_mat(fine_obj->myhexahedrons[fine_obj->myhexahedrons.size()-1],PoissonRatio,YoungModulus1);
		}
	    }	
	}
    }
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  i=1;
  stringstream ss_fine; string sub_fine; ss_fine<<i; ss_fine>>sub_fine;
  string path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".vtk";
  myio.saveAsVTK(fine_obj,path_name_fine);
  string mat_path_name_fine=out_dir+"/"+object_name+".sub"+sub_fine+".mat";
  myio.saveMatPara(fine_obj,mat_path_name_fine,1);// 1表示save material_para array ; 0表示save stiffness_tensor;
  delete fine_obj;
  return 0;
}

int object_creator::bind_mat(hexahedron &hexahedron_here,const double &PoissonRatio,const double &YoungModulus)
{
  size_t i,j,k;
  double coefficient1,coefficient2;
  coefficient1=YoungModulus/(1+PoissonRatio)*0.5;
  coefficient2=YoungModulus*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);
  for(i=0;i<2;++i)
    {
      for(j=0;j<2;++j)
	{
	  for(k=0;k<2;++k)
	    {
	      hexahedron_here.material_para[i][j][k][0]=coefficient1;
	      hexahedron_here.material_para[i][j][k][1]=coefficient2;
	    }
	}
    }
  return 0;
}
