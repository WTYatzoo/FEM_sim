#include "CtoFer.h"
#include "vertex.h"
#include "hexahedron.h"
#include "io.h"
#include "util.h"
using namespace std;
using namespace Eigen;
const static double help[8][3]= {
  {
    -1,-1,-1
  },
  {
    -1,-1,1
  },
  {
    -1,1,-1
  },
  {
    -1,1,1
  },
  {
    1,-1,-1
  },
  {
    1,-1,1
  },
  {
    1,1,-1
  },
  {
    1,1,1
  }
};

const static double quadrature[1][3]={
  {
    0, 0, 0
  }
};

CtoFer::CtoFer(boost::property_tree::ptree &para_tree)
{
  io myio=io();
  string input_fine_object=para_tree.get<string>("input_fine_object.value");
  string input_coarsen_deform_object=para_tree.get<string>("input_coarsen_deform_object.value");
  string input_coarsen_rest_object=para_tree.get<string>("input_coarsen_rest_object.value");
  string input_harmonic_def=para_tree.get<string>("input_harmonic_def.value"); // harmonic displacement 的下标从1开始，不是从0开始
  string object_name=para_tree.get<string>("object_name.value");
  string out_dir=para_tree.get<string>("out_dir.value");

  int level=para_tree.get<int>("level.value"); //for edge: 1 to 2 ; 1 to 4 or 1 to 8
  
  size_t i,j,k;
  size_t a,b,c;
  size_t x,y,z;

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
  
  object* coarsen_obj=new object();
  myio.getVertexAndHex(coarsen_obj->myvertexs,coarsen_obj->myhexahedrons,input_coarsen_rest_object);
  coarsen_obj->num_hexahedrons=coarsen_obj->myhexahedrons.size();
  coarsen_obj->num_vertex=coarsen_obj->myvertexs.size();

  vector<vertex > myvertexs_deform;
  vector<hexahedron > myhexahedrons_deform;
  myio.getVertexAndHex(myvertexs_deform,myhexahedrons_deform,input_coarsen_deform_object);
  for(i=0;i<coarsen_obj->num_vertex;i++)
    {
      coarsen_obj->myvertexs[i].location=coarsen_obj->myvertexs[i].location_maybe=myvertexs_deform[i].location;
    }

  double dmetric_coarsen=myio.getDmetric(input_coarsen_rest_object);
  double dmetric_fine=dmetric_coarsen/(double)level;
  printf("dmetric_coarsen:: %lf\ndmetric_fine:: %lf\n",dmetric_coarsen,dmetric_fine);
  
  size_t num_hex=coarsen_obj->myhexahedrons.size();
  size_t num_fine_in_coarsen=pow(level,3);
  printf("coarsen num_hex:: %d\n",num_hex);

  vector<vertex > myvertexs_fine;
  vector<hexahedron > myhexahedrons_fine;

  size_t ct_vertex=0;
  size_t ct_hex=0;
  size_t index_vertex_new[2][2][2];

  double shapeFuncGrad[2][2][2][3];
  calShapeFuncGrad(shapeFuncGrad);


  for(i=0;i<num_hex;++i) //coarsen hex index
    {
      MatrixXd R(3,3);  MatrixXd S(3,3);
      calRandS(dmetric_coarsen,R,S,coarsen_obj->myvertexs,coarsen_obj->myhexahedrons[i],shapeFuncGrad);
      MatrixXd bc_ori(3,1); MatrixXd bc(3,1);
      calBaryCenter(bc_ori,bc,coarsen_obj->myvertexs,coarsen_obj->myhexahedrons[i]);
      int ct=0;
      MatrixXd H(3,6); MatrixXd H_bc(3,6); VectorXd S_array(6); MatrixXd col_H_bc(3,1); MatrixXd col_H(3,1);

      for(a=0;a<3;a++)
	{
	  for(b=a;b<3;b++)
	    {
	      //   calHbc(col_H_bc,myvertexs[ct],myhexahedrons[ct],i,num_fine_in_coarsen,bc_ori);
	      calHbc(col_H_bc,myvertexs[ct],coarsen_obj->myhexahedrons[i],bc_ori);
	      if(a==b)
		{
		  H_bc.col(ct)=col_H_bc;
		}
	      else if(a!=b)
		{
		  H_bc.col(ct)=col_H_bc;
		}
	      ct++;
	    }
	}


      MatrixXd A(24,6); VectorXd B(24,1);
      for(a=0;a<8;a++)
	{
	  x=a/4; y=(a%4)/2; z=(a%4)%2;
	  size_t index_vertex_here=coarsen_obj->myhexahedrons[i].index_vertex[x][y][z];
	  ct=0;
	  for(x=0;x<3;x++)
	    {
	      for(y=x;y<3;y++)
		{
		  calH(col_H,myvertexs[ct],fine_obj,index_vertex_here);
		  if(x==y)
		    {
		      H.col(ct)=col_H;
		    }
		  else if(x!=y)
		    {
		      H.col(ct)=col_H;
		    }		      
		  ct++;
		}
	    }
	  for(b=a*3;b<a*3+3;b++)
	    {
	      A.row(b)=(H-H_bc).row(b%3);
	    }

	  MatrixXd loc_now(3,1); MatrixXd loc_ori(3,1);
	  for(b=0;b<3;b++)
	    {
	      loc_ori(b,0)=coarsen_obj->myvertexs[index_vertex_here].location_original(b);
	      loc_now(b,0)=coarsen_obj->myvertexs[index_vertex_here].location(b);
	    }
	  VectorXd B_here=R.transpose()*(loc_now-bc)-loc_ori+bc_ori;
	  for(b=a*3;b<a*3+3;b++)
	    {
	      B(b)=B_here(b%3);
	    }
	}

      MatrixXd ATA=A.transpose()*A;
      VectorXd ATB=A.transpose()*B;
      S_array=ATA.llt().solve(ATB);
      cout<<"S"<<endl<<S_array<<endl;
      
      for(j=i*num_fine_in_coarsen;j<i*num_fine_in_coarsen+num_fine_in_coarsen;j++)
	{
	  for(a=0;a<2;a++)
	    {
	      for(b=0;b<2;b++)
		{
		  for(c=0;c<2;c++)
		    {
		      size_t index_vertex_here=fine_obj->myhexahedrons[j].index_vertex[a][b][c];
		      ct=0;
		      for(x=0;x<3;x++)
			{
			  for(y=x;y<3;y++)
			    {
			      calH(col_H,myvertexs[ct],fine_obj,index_vertex_here);
			      if(x==y)
				{
				  H.col(ct)=col_H;
				}
			      else if(x!=y)
				{
				  H.col(ct)=col_H;
				}
			      
			      ct++;
			    }
			}

		      MatrixXd loc_ori(3,1);
		      for(x=0;x<3;x++)
			{
			  loc_ori(x,0)=fine_obj->myvertexs[index_vertex_here].location(x);
			}
		      MatrixXd x_new=R*(loc_ori-bc_ori+(H-H_bc)*S_array)+bc;
		      myvector loc_new;
		      for(x=0;x<3;x++)
			{
			  loc_new(x)=x_new(x,0);
			}
		      vertex vt_new=vertex(loc_new);
		      myvertexs_fine.push_back(vt_new);
		      index_vertex_new[a][b][c]=ct_vertex;
		      ct_vertex++;		      
		    }
		}
	    }
	  hexahedron hex_new=hexahedron(index_vertex_new);
	  myhexahedrons_fine.push_back(hex_new);
	  ct_hex++;
	}
    }

  /*
  for(i=0;i<num_hex;++i) //coarsen hex index
    {
      MatrixXd R(3,3);  MatrixXd S(3,3);
      calRandS(dmetric_coarsen,R,S,coarsen_obj->myvertexs,coarsen_obj->myhexahedrons[i],shapeFuncGrad);
      cout<<R<<endl; cout<<S<<endl;

      
      MatrixXd bc_ori(3,1); MatrixXd bc(3,1);
      calBaryCenter(bc_ori,bc,coarsen_obj->myvertexs,coarsen_obj->myhexahedrons[i]);
      int ct=0;
      MatrixXd H(3,6); MatrixXd H_bc(3,6); MatrixXd S_array(6,1); MatrixXd col_H_bc(3,1); MatrixXd col_H(3,1);
      for(a=0;a<3;a++)
	{
	  for(b=a;b<3;b++)
	    {
	      calHbc(col_H_bc,myvertexs[ct],myhexahedrons[ct],i,num_fine_in_coarsen,bc_ori);

	      if(a==b)
		{
		  H_bc.col(ct)=col_H_bc;
		}
	      else if(a!=b)
		{
		  H_bc.col(ct)=2*col_H_bc;
		}

	      if(a==b)
		{
		  S_array(ct,0)=S(a,b);
		}
	      else if(a!=b)
		{
		  S_array(ct,0)=2*S(a,b);
		}
	      ct++;
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
		      size_t index_vertex_here=fine_obj->myhexahedrons[j].index_vertex[a][b][c];
		      ct=0;
		      for(x=0;x<3;x++)
			{
			  for(y=x;y<3;y++)
			    {
			      calH(col_H,myvertexs[ct],fine_obj,index_vertex_here);
			      if(x==y)
				{
				  H.col(ct)=col_H;
				}
			      else if(x!=y)
				{
				  H.col(ct)=2*col_H;
				}
			      
			      ct++;
			    }
			}

		      MatrixXd loc_ori(3,1);
		      for(x=0;x<3;x++)
			{
			  loc_ori(x,0)=fine_obj->myvertexs[index_vertex_here].location(x);
			}
		      MatrixXd x_new=R*(loc_ori-bc_ori+(H-H_bc)*S_array)+bc;
		      myvector loc_new;
		      for(x=0;x<3;x++)
			{
			  loc_new(x)=x_new(x,0);
			}
		      vertex vt_new=vertex(loc_new);
		      myvertexs_fine.push_back(vt_new);
		      index_vertex_new[a][b][c]=ct_vertex;
		      ct_vertex++;		      
		    }
		}
	    }
	  hexahedron hex_new=hexahedron(index_vertex_new);
	  myhexahedrons_fine.push_back(hex_new);
	  ct_hex++;
	}
	}
  */
  /*
   for(i=0;i<num_hex;++i) //coarsen hex index
    {
      MatrixXd R(3,3);  MatrixXd S(3,3);
      calRandS(dmetric_coarsen,R,S,coarsen_obj->myvertexs,coarsen_obj->myhexahedrons[i],shapeFuncGrad);
      cout<<R<<endl; cout<<S<<endl;

      
      MatrixXd center_ori(3,1); MatrixXd center(3,1);

      for(a=0;a<3;a++)
	{
	  center_ori(a,0)=coarsen_obj->myvertexs[coarsen_obj->myhexahedrons[i].index_vertex[0][0][0]].location_original(a);
	  center(a,0)=coarsen_obj->myvertexs[coarsen_obj->myhexahedrons[i].index_vertex[0][0][0]].location(a);
	}
      int ct=0;
      MatrixXd H(3,6); MatrixXd H_c(3,6); MatrixXd S_array(6,1); MatrixXd col_H_c(3,1); MatrixXd col_H(3,1);
      for(a=0;a<3;a++)
	{
	  for(b=a;b<3;b++)
	    {
	      calH(col_H_c,myvertexs[ct],fine_obj,coarsen_obj->myhexahedrons[i].index_vertex[0][0][0]);
	      H_c.col(ct)=col_H_c;
	      S_array(ct,0)=S(a,b);
	      ct++;
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
		      size_t index_vertex_here=fine_obj->myhexahedrons[j].index_vertex[a][b][c];
		      ct=0;
		      for(x=0;x<3;x++)
			{
			  for(y=x;y<3;y++)
			    {
			      calH(col_H,myvertexs[ct],fine_obj,index_vertex_here);
			      H.col(ct)=col_H;
			      ct++;
			    }
			}

		      MatrixXd loc_ori(3,1);
		      for(x=0;x<3;x++)
			{
			  loc_ori(x,0)=fine_obj->myvertexs[index_vertex_here].location(x);
			}
		      MatrixXd x_new=R*(loc_ori-center_ori+2*(H-H_c)*S_array)+center;
		      myvector loc_new;
		      for(x=0;x<3;x++)
			{
			  loc_new(x)=x_new(x,0);
			}
		      vertex vt_new=vertex(loc_new);
		      myvertexs_fine.push_back(vt_new);
		      index_vertex_new[a][b][c]=ct_vertex;
		      ct_vertex++;		      
		    }
		}
	    }
	  hexahedron hex_new=hexahedron(index_vertex_new);
	  myhexahedrons_fine.push_back(hex_new);
	  ct_hex++;
	}
    }
  */
  string fine_obj_from_coarsen=input_coarsen_deform_object+".refine.vtk";
  //  string fine_obj_from_coarsen=out_dir+"/"+object_name+".sub1.new.vtk"; 
  myio.saveArrayAsVTK(myvertexs_fine,myhexahedrons_fine,fine_obj_from_coarsen);  
  delete coarsen_obj;
  delete fine_obj;
}

int CtoFer::calH(MatrixXd &col_H,vector<vertex > myvertexs,object *fine_obj,size_t index_vertex_here)
{
  size_t i;
  for(i=0;i<3;i++)
    {
      col_H(i,0)=myvertexs[index_vertex_here].location(i)-fine_obj->myvertexs[index_vertex_here].location_original(i);
    }
  return 0;
}

int CtoFer::calHbc(MatrixXd &col_H_bc,vector<vertex > &myvertexs,const vector<hexahedron> &myhexahedrons,size_t i,size_t num_fine_in_coarsen,MatrixXd &bc_ori)
 {
   size_t j;
   size_t a,b,c;
   map<int,int > mp; mp.clear();

   size_t ct=0;
   myvector help=myvector(0,0,0);
   for(j=i*num_fine_in_coarsen;j<i*num_fine_in_coarsen+num_fine_in_coarsen;j++)
     {
       for(a=0;a<2;a++)
	 {
	   for(b=0;b<2;b++)
	     {
	       for(c=0;c<2;c++)
		 {
		   int index_vertex_here=myhexahedrons[j].index_vertex[a][b][c];
		   if(mp.find(index_vertex_here)==mp.end())
		     {
		       help+=myvertexs[index_vertex_here].location;
		       ct++;
		       mp[index_vertex_here]=1;
		     }
		 }
	     }
	 }
     }

   cout<<"vertex in fine of coarsen:"<<ct<<endl;
   help/=(double)ct;
   for(j=0;j<3;j++)
     {
       col_H_bc(j,0)=help(j)-bc_ori(j,0);
     }
   return 0;
 }

int CtoFer::calHbc(MatrixXd &col_H_bc,std::vector<vertex > &myvertexs,const hexahedron &hexahedron_here,MatrixXd &bc_ori)
{
  size_t a,b,c,d;
  col_H_bc.fill(0);
  for(a=0;a<2;a++)
    {
      for(b=0;b<2;b++)
	{
	  for(c=0;c<2;c++)
	    {
	      size_t index_vertex_here=hexahedron_here.index_vertex[a][b][c];
	      for(d=0;d<3;d++)
		{
		  col_H_bc(d,0)+=myvertexs[index_vertex_here].location(d);
		}
	    }	  
	}
    }
  col_H_bc*=0.125;
  col_H_bc-=bc_ori;
  return 0;
}
int CtoFer::calBaryCenter(MatrixXd &bc_ori,MatrixXd &bc, std::vector<vertex > &myvertexs,const hexahedron &hexahedron_here)
{
  size_t a,b,c,d;
  bc_ori.fill(0); bc.fill(0);
  for(a=0;a<2;a++)
    {
      for(b=0;b<2;b++)
	{
	  for(c=0;c<2;c++)
	    {
	      for(d=0;d<3;d++)
		{
		  bc_ori(d,0)+=myvertexs[hexahedron_here.index_vertex[a][b][c]].location_original(d);
		  bc(d,0)+=myvertexs[hexahedron_here.index_vertex[a][b][c]].location(d);
		}	     
	    }
	}
    }
  bc_ori*=0.125; bc*=0.125;
  return 0;
}
int CtoFer::calRandS(const double &dmetric,Eigen::MatrixXd &R,Eigen::MatrixXd &S, std::vector< vertex > &myvertexs,const hexahedron &hexahedron_here,const double (&shapeFuncGrad)[2][2][2][3])
  
{
  size_t a,b,c;
  size_t i,j,k;
 
  MatrixXd F(3,3);  
  size_t row,col;
  
  double help=2.0/dmetric;       
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
		      value_now+=myvertexs[hexahedron_here.index_vertex[i][j][k]].location(row)*shapeFuncGrad[i][j][k][col];
		    }
			      
		}
	    }
	  F(row,col)=value_now;
	}
    }
  F*=help; //deformation gradient

  util myutil=util();
  myutil.extractRotation(F,R);
  // JacobiSVD<MatrixXd> svd(F, ComputeThinU | ComputeThinV);
  // R=svd.matrixU()*(svd.matrixV().adjoint());
  S=R.transpose()*F;
  return 0;
}

int CtoFer::calShapeFuncGrad(double (&shapeFuncGrad)[2][2][2][3])
{
  size_t i,j,k;
  size_t whichShapeFunc;
  size_t whichQuadrature=0;
  for(i=0;i<2;++i)
    {
      for(j=0;j<2;++j)
	{
	  for(k=0;k<2;++k)
	    {
	      whichShapeFunc=4*i+2*j+k;	      
	      shapeFuncGrad[i][j][k][0]=help[whichShapeFunc][0]*(1+help[whichShapeFunc][1]*quadrature[whichQuadrature][1])*(1+help[whichShapeFunc][2]*quadrature[whichQuadrature][2])*0.125;
	      shapeFuncGrad[i][j][k][1]=help[whichShapeFunc][1]*(1+help[whichShapeFunc][0]*quadrature[whichQuadrature][0])*(1+help[whichShapeFunc][2]*quadrature[whichQuadrature][2])*0.125;
	      shapeFuncGrad[i][j][k][2]=help[whichShapeFunc][2]*(1+help[whichShapeFunc][1]*quadrature[whichQuadrature][1])*(1+help[whichShapeFunc][0]*quadrature[whichQuadrature][0])*0.125;
	    }
	}
    }
  return 0;
}
