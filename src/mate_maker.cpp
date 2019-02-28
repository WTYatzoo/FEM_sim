#include "mate_maker.h"
#include "io.h"
#include "object.h"
using namespace std;

mate_maker::mate_maker(boost::property_tree::ptree &para_tree)
{
  string file_vtk_fine=para_tree.get<string >("file_vtk_fine.value");
  string mat_out_fine=para_tree.get<string >("mat_out_fine.value");
  string mat_out_fine_txt=para_tree.get<string>("mat_out_fine_txt.value");
  string vtk_out_fine=para_tree.get<string >("vtk_out_fine.value");
  double y[2]; y[0]=para_tree.get<double >("y1.value");  y[1]=para_tree.get<double >("y2.value");
  double p=para_tree.get<double >("p.value");
  int axis=para_tree.get<int >("axis.value"); // 0:x 1:y 2:z 哪一维是重复的

  io myio=io();
  double dmetric=myio.getDmetric(file_vtk_fine);
  printf("dmetric:: %lf\n",dmetric);
  object* fine_obj=new object();
  myio.getVertexAndHex(fine_obj->myvertexs,fine_obj->myhexahedrons,file_vtk_fine);
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  fine_obj->num_vertex=fine_obj->myvertexs.size();

  map<pair<int, int >,int > mp_go; mp_go.clear();
  map<int,pair<int, int > > mp_come; mp_come.clear();
  map<int,int > mp_y; mp_y.clear();
  
  int num_hexahedrons=fine_obj->num_hexahedrons;
  vector<int > index_hex[num_hexahedrons];

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
  myio.saveMatPara(fine_obj,mat_out_fine,1);
  myio.saveAsVTKwithPara(fine_obj,vtk_out_fine);
  myio.saveMatParaAsTXT(fine_obj,mat_out_fine_txt);
  
  delete fine_obj;
}

int mate_maker::cal_mate(hexahedron &hex_now,double y,double p)
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
