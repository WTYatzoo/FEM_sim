#include "mate_factory.h"
#include "io.h"
#include "object.h"
using namespace std;

vector<int > index_hex[700000];

mate_factory::mate_factory(boost::property_tree::ptree &para_tree)
{
  string file_vtk_fine=para_tree.get<string >("file_vtk_fine.value");
  string file_out_fine=para_tree.get<string >("file_out_fine.value");
  double y_homo=para_tree.get<double >("y_homo.value");
  double y[2]; y[0]=para_tree.get<double >("y1.value");  y[1]=para_tree.get<double >("y2.value");
  double p=para_tree.get<double >("p.value");
  string object_name=para_tree.get<string >("object_name.value");

  double N=para_tree.get<double >("N.value");
  double C=para_tree.get<double >("C.value");

  io myio=io();
  dmetric=myio.getDmetric(file_vtk_fine);
  printf("dmetric:: %lf\n",dmetric);

  object* fine_obj=new object();
  myio.getVertexAndHex(fine_obj->myvertexs,fine_obj->myhexahedrons,file_vtk_fine);
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  
  size_t i,j;
  int ct_mate=0;
  for(i=0;i<7;i++)
    {
      if(i==0)
	{
	  mate_A(fine_obj,y_homo,p);
	  stringstream ss_ct_mate; string str_ct_mate; ss_ct_mate<<ct_mate; ss_ct_mate>>str_ct_mate;
	  string vtk_out_fine=file_out_fine+"/"+object_name+"_"+str_ct_mate+".vtk";
	  string mat_out_fine_txt=file_out_fine+"/"+object_name+"_"+str_ct_mate+".txt";
	  myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
	  myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
	  ct_mate++;
	}
      else if(i==1)
	{
	  for(j=0;j<3;j++)
	    {
	      mate_B(fine_obj,y,p,j);
	      stringstream ss_ct_mate; string str_ct_mate; ss_ct_mate<<ct_mate; ss_ct_mate>>str_ct_mate;
	      string vtk_out_fine=file_out_fine+"/"+object_name+"_"+str_ct_mate+".vtk";
	      string mat_out_fine_txt=file_out_fine+"/"+object_name+"_"+str_ct_mate+".txt";
	      myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
	      myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
	      ct_mate++;
	    }
	}
      else if(i==2)
	{
	  for(j=0;j<3;j++)
	    {
	      mate_C(fine_obj,y,p,j);
	      stringstream ss_ct_mate; string str_ct_mate; ss_ct_mate<<ct_mate; ss_ct_mate>>str_ct_mate;
	      string vtk_out_fine=file_out_fine+"/"+object_name+"_"+str_ct_mate+".vtk";
	      string mat_out_fine_txt=file_out_fine+"/"+object_name+"_"+str_ct_mate+".txt";
	      myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
	      myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
	      ct_mate++;
	    }
	}
      else if(i==3)
	{
	  mate_D(fine_obj,y,p);
	  stringstream ss_ct_mate; string str_ct_mate; ss_ct_mate<<ct_mate; ss_ct_mate>>str_ct_mate;
	  string vtk_out_fine=file_out_fine+"/"+object_name+"_"+str_ct_mate+".vtk";
	  string mat_out_fine_txt=file_out_fine+"/"+object_name+"_"+str_ct_mate+".txt";
	  myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
	  myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
	  ct_mate++;
	}
      else if(i==4)
	{
	  for(j=0;j<3;j++)
	    {
	      mate_E(fine_obj,y,p,j);
	      stringstream ss_ct_mate; string str_ct_mate; ss_ct_mate<<ct_mate; ss_ct_mate>>str_ct_mate;
	      string vtk_out_fine=file_out_fine+"/"+object_name+"_"+str_ct_mate+".vtk";
	      string mat_out_fine_txt=file_out_fine+"/"+object_name+"_"+str_ct_mate+".txt";
	      myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
	      myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
	      ct_mate++;
	    }  
	  
	}
      else if(i==5)
	{
	  mate_F(fine_obj,y,p);
	  stringstream ss_ct_mate; string str_ct_mate; ss_ct_mate<<ct_mate; ss_ct_mate>>str_ct_mate;
	  string vtk_out_fine=file_out_fine+"/"+object_name+"_"+str_ct_mate+".vtk";
	  string mat_out_fine_txt=file_out_fine+"/"+object_name+"_"+str_ct_mate+".txt";
	  myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
	  myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
	  ct_mate++;
	}
      else if(i==6)
	{
	  mate_G(fine_obj,y,p,N,C);
	  stringstream ss_ct_mate; string str_ct_mate; ss_ct_mate<<ct_mate; ss_ct_mate>>str_ct_mate;
	  string vtk_out_fine=file_out_fine+"/"+object_name+"_"+str_ct_mate+".vtk";
	  string mat_out_fine_txt=file_out_fine+"/"+object_name+"_"+str_ct_mate+".txt";
	  myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
	  myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
	  ct_mate++;
	}
    }
  delete fine_obj;
}

int mate_factory::mate_A(object* fine_obj,double y_homo,double p)
{
  size_t num_hexahedrons=fine_obj->num_hexahedrons;
  for(size_t i=0;i<num_hexahedrons;i++)
    {
      cal_mate(fine_obj->myhexahedrons[i],y_homo,p);
    }
  return 0;
}

int mate_factory::mate_B(object* fine_obj,double (&y)[2],double p,int axis)
{
  if(y[0]<y[1])
    {
      double help=y[0]; y[0]=y[1]; y[1]=help; //保证y[0]是大的
    }  
  int num_hexahedrons=fine_obj->num_hexahedrons;
  size_t i,j,k;
  size_t use[1]; use[0]=axis;
  double loc_min[1];

  for(j=0;j<1;j++)
    {
      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[0].index_vertex[0][0][0]].location(use[j]);
    }
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<1;j++)
	{
	  if(fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j])<loc_min[j])
	    {
	      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j]);
	    }
	}
    }

  int index[1];
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<1;j++)
	{
	  index[j]=(int)((fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j])-loc_min[j])/dmetric+0.5);
	}
      cal_mate(fine_obj->myhexahedrons[i],y[index[j]%2],p);
    }
  return 0;
}

int mate_factory::mate_C(object* fine_obj,double (&y)[2],double p,int axis)
{
  if(y[0]<y[1])
    {
      double help=y[0]; y[0]=y[1]; y[1]=help; //保证y[0]是大的
    }  
  int num_hexahedrons=fine_obj->num_hexahedrons;
  size_t i,j,k;
  size_t use[1]; use[0]=axis;
  double loc_min[1];

  for(j=0;j<1;j++)
    {
      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[0].index_vertex[0][0][0]].location(use[j]);
    }
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<1;j++)
	{
	  if(fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j])<loc_min[j])
	    {
	      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j]);
	    }
	}
    }

  int index[1];
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<1;j++)
	{
	  index[j]=(int)((fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j])-loc_min[j])/dmetric+0.5);
	}
      if(index[j]%4==1||index[j]%4==2)
	{
	  cal_mate(fine_obj->myhexahedrons[i],y[1],p);
	}
      else if(index[j]%4==0||index[j]%4==3)
	{
	  cal_mate(fine_obj->myhexahedrons[i],y[0],p);
	}
      
    }
  return 0;
}

int mate_factory::mate_D(object* fine_obj,double (&y)[2],double p)
{
  int num_hexahedrons=fine_obj->num_hexahedrons;
  size_t i,j,k;
  double loc_min[3];

  for(j=0;j<1;j++)
    {
      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[0].index_vertex[0][0][0]].location(j);
    }
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<3;j++)
	{
	  if(fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(j)<loc_min[j])
	    {
	      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(j);
	    }
	}
    }

  int index[3];
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<3;j++)
	{
	  index[j]=(int)((fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(j)-loc_min[j])/dmetric+0.5);
	}
      int index_mod[3];
      for(j=0;j<3;j++)
	{
	  index_mod[j]=index[j]%2;
	}
      int help=index_mod[0]^index_mod[1]^index_mod[2];
      cal_mate(fine_obj->myhexahedrons[i],y[help],p);      
    }
  return 0;
}


int mate_factory::mate_E(object* fine_obj,double (&y)[2],double p,int axis)
{
  map<pair<int, int >,int > mp_go; mp_go.clear();
  map<int,pair<int, int > > mp_come; mp_come.clear();
  map<int,int > mp_y; mp_y.clear();
  
  int num_hexahedrons=fine_obj->num_hexahedrons;
  printf("fine_obj num_hex: %d\n",num_hexahedrons);
  
  size_t i,j,k;
  size_t a,b;
  size_t use[2];
  if(axis==0)
    {
      use[0]=1; use[1]=2; 
    }
  else if(axis==1)
    {
      use[0]=0;  use[1]=2;
    }
  else if(axis==2)
    {
      use[0]=0; use[1]=1;
    }

  for(i=0;i<num_hexahedrons;i++)
    {
      while(!index_hex[i].empty())
	{
	  index_hex[i].pop_back();
	}
    }
  
  double loc_min[2];

  for(j=0;j<2;j++)
    {
      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[0].index_vertex[0][0][0]].location(use[j]);
    }
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<2;j++)
	{
	  if(fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j])<loc_min[j])
	    {
	      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j]);
	    }
	}
    }

  int ct=0;
  int index[2];
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<2;j++)
	{
	  index[j]=(int)((fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(use[j])-loc_min[j])/dmetric+0.5);
	}
      if(mp_go.find(make_pair(index[0],index[1]))==mp_go.end())
	{
	  mp_go[make_pair(index[0],index[1])]=ct;
	  mp_come[ct]=make_pair(index[0],index[1]);
	  index_hex[ct].push_back(i);
	  ct++;
	}
      else if(mp_go.find(make_pair(index[0],index[1]))!=mp_go.end())
	{
	  int ct_now=mp_go[make_pair(index[0],index[1])];
	  index_hex[ct_now].push_back(i);
	}      
    }

  queue<int > que;
  while(!que.empty())
    {
      que.pop();
    }
  que.push(0);
  mp_y[0]=1;

  
  int loc_now[2];
  loc_now[0]=mp_come[0].first; loc_now[1]=mp_come[0].second;
  for(a=0;a<=1;a++)
    {
      int loc_here[2];
      loc_here[0]=loc_now[0]+a*2-1;
      loc_here[1]=loc_now[1];

      if(mp_go.find(make_pair(loc_here[0],loc_here[1]))!=mp_go.end())
	{
	  int which=mp_go[make_pair(loc_here[0],loc_here[1])];
	  mp_y[which]=mp_y[0]^a;
	  que.push(which);	   
	}
    }
  for(b=0;b<=1;b++)
    {
      int loc_here[2];
      loc_here[0]=loc_now[0];
      loc_here[1]=loc_now[1]+b*2-1;
      if(mp_go.find(make_pair(loc_here[0],loc_here[1]))!=mp_go.end())
	{
	  int which=mp_go[make_pair(loc_here[0],loc_here[1])];
	  mp_y[which]=mp_y[0]^b;
	  que.push(which);
	}
    }
  
  while(!que.empty())
    {
      int siz=que.size();
      while(siz--)
	{
	  int now=que.front();
	  que.pop();
	  loc_now[0]=mp_come[now].first;
	  loc_now[1]=mp_come[now].second;
	  bool inverse;
	  int find=0;
	  for(a=0;a<=1;a++)
	    {
	      int loc_here[2];
	      loc_here[0]=loc_now[0]+2*a-1;
	      loc_here[1]=loc_now[1];
	      if(mp_go.find(make_pair(loc_here[0],loc_here[1]))!=mp_go.end())
		{
		  int which=mp_go[make_pair(loc_here[0],loc_here[1])];
		  if(mp_y.find(which)!=mp_y.end()&&find==0)
		    {
		      inverse=a^mp_y[now]^mp_y[which];
		      find=1;
		    }		 
		}
	    }
	  for(b=0;b<=1;b+=1)
	    {
	      int loc_here[2];
	      loc_here[0]=loc_now[0];
	      loc_here[1]=loc_now[1]+2*b-1;
	      if(mp_go.find(make_pair(loc_here[0],loc_here[1]))!=mp_go.end())
		{
		  int which=mp_go[make_pair(loc_here[0],loc_here[1])];
		  if(mp_y.find(which)!=mp_y.end()&&find==0)
		    {
		      inverse=b^mp_y[now]^mp_y[which];
		      find=1;
		    }		 
		}
	    }
	  for(a=0;a<=1;a+=1)
	    {
	      int loc_here[2];
	      loc_here[0]=loc_now[0]+a*2-1;
	      loc_here[1]=loc_now[1];
	      if(mp_go.find(make_pair(loc_here[0],loc_here[1]))!=mp_go.end())
		{
		  int which=mp_go[make_pair(loc_here[0],loc_here[1])];
		  if(mp_y.find(which)==mp_y.end())
		    {
		      mp_y[which]=mp_y[now]^a^inverse;
		      que.push(which);
		    }		 
		}
	    }
	  for(b=0;b<=1;b++)
	    {
	      int loc_here[2];
	      loc_here[0]=loc_now[0];
	      loc_here[1]=loc_now[1]+b*2-1;
	      if(mp_go.find(make_pair(loc_here[0],loc_here[1]))!=mp_go.end())
		{
		  int which=mp_go[make_pair(loc_here[0],loc_here[1])];
		  if(mp_y.find(which)==mp_y.end())
		    {
		      mp_y[which]=mp_y[now]^b^inverse;
		      que.push(which);		      
		    }		 
		}
	    }	  
	}      
    }
  
  for(i=0;i<ct;i++)
    {
      vector<int>::iterator it;
      for(it=index_hex[i].begin();it!=index_hex[i].end();it++)
	{
	  size_t index_hex_now=(*it);
	  if(mp_y.find(i)==mp_y.end())
	    {
	      cal_mate(fine_obj->myhexahedrons[index_hex_now],0.0,p);
	    }
	  else
	    {
	      int index_y_now=mp_y[i];
	      cal_mate(fine_obj->myhexahedrons[index_hex_now],y[index_y_now],p);
	    }
	}
	}
  return 0;
}

int mate_factory::mate_F(object* fine_obj,double (&y)[2],double p)
{
  if(y[0]<y[1])
    {
      double help=y[0]; y[0]=y[1]; y[1]=help; //保证y[0]是大的
    }  
  int num_hexahedrons=fine_obj->num_hexahedrons;
  size_t i,j,k;  
  double loc_min[3];

  for(j=0;j<2;j++)
    {
      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[0].index_vertex[0][0][0]].location(j);
    }
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<3;j++)
	{
	  if(fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(j)<loc_min[j])
	    {
	      loc_min[j]=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(j);
	    }
	}
    }

  int index[3];
  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      for(j=0;j<3;j++)
	{
	  index[j]=(int)((fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[0][0][0]].location(j)-loc_min[j])/dmetric+0.5);
	}
      int index_mod[3];  bool check[3];
      for(j=0;j<3;j++)
	{
	  index_mod[j]=index[j]%4;
	  check[j]=(index_mod[j]==1||index_mod[j]==2);
	}
      if(check[0]&&check[1]&&check[2])
	{
	  cal_mate(fine_obj->myhexahedrons[i],y[1],p);
	}
      else
	{
	  cal_mate(fine_obj->myhexahedrons[i],y[0],p);
	}
    }
  return 0;
}

int mate_factory::mate_G(object* fine_obj,double (&y)[2],double p,double N,double C)
// level set function is sin(N*pi*x)*cos(N*pi*y) + sin(N*pi*z)*cos(N*pi*x) + sin(N*pi*y)*cos(N*pi*z) - C
// where negative is hard and positive is soft
{
  size_t i,j,k;
  size_t a,b,c;
  if(y[0]<y[1])
    {
      double help=y[0]; y[0]=y[1]; y[1]=help; //保证y[0]是大的
    }

  myvector loc_now;
  for(i=0;i<fine_obj->num_vertex;i++)
    {
      loc_now=fine_obj->myvertexs[i].location;
      fine_obj->myvertexs[i].scalar_field_value=sin(N*pi*loc_now(0))*cos(N*pi*loc_now(1)) + sin(N*pi*loc_now(2))*cos(N*pi*loc_now(0)) + sin(N*pi*loc_now(1))*cos(N*pi*loc_now(2)) - C;
    }

  for(i=0;i<fine_obj->num_hexahedrons;i++)
    {
      double sum=0;
      for(a=0;a<2;a++)
	{
	  for(b=0;b<2;b++)
	    {
	      for(c=0;c<2;c++)
		{
		  sum+=fine_obj->myvertexs[fine_obj->myhexahedrons[i].index_vertex[a][b][c]].scalar_field_value;
		}
	    }
	}
      if(sum<1e-16)
	{
	  cal_mate(fine_obj->myhexahedrons[i],y[0],p);
	}
      else
	{
	  cal_mate(fine_obj->myhexahedrons[i],y[1],p);
	}
    }
  return 0;
}
int mate_factory::cal_mate(hexahedron &hex_now,double y,double p)
{
  size_t i,j,k;
  double c1=y/(1+p)*0.5;
  double c2=y*p/(1+p)/(1-2*p);
  for(i=0;i<2;i++)
    {
      for(j=0;j<2;j++)
	{
	  for(k=0;k<2;k++)
	    {	      
	      hex_now.material_para[i][j][k][0]=c1;
	      hex_now.material_para[i][j][k][1]=c2;
	    }
	}
    }
  return 0;
}
