#include "myvector.h"
#include "error_measurer.h"
#include "io.h"
using namespace std;
error_measurer::error_measurer(boost::property_tree::ptree &para_tree)
{
  size_t i,j;
  string fine_deform_mesh=para_tree.get<string >("fine_deform_mesh.value");
  string coarsen_deform_mesh=para_tree.get<string >("coarsen_deform_mesh.value");
  string max_error_txt=para_tree.get<string >("max_error_txt.value");
  string average_error_txt=para_tree.get<string >("average_error_txt.value");
  
  io myio=io();
  object* fine_obj=new object();
  myio.getVertexAndHex(fine_obj->myvertexs,fine_obj->myhexahedrons,fine_deform_mesh);
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  fine_obj->num_vertex=fine_obj->myvertexs.size();

  object* coarsen_obj=new object();
  myio.getVertexAndHex(coarsen_obj->myvertexs,coarsen_obj->myhexahedrons,coarsen_deform_mesh);
  coarsen_obj->num_hexahedrons=coarsen_obj->myhexahedrons.size();
  coarsen_obj->num_vertex=coarsen_obj->myvertexs.size();

  size_t num_vertex_coarsen=coarsen_obj->num_vertex;
  size_t num_vertex_fine=fine_obj->num_vertex;

  myvector loc_max=fine_obj->myvertexs[0].location;
  myvector loc_min=fine_obj->myvertexs[0].location;

  for(i=0;i<num_vertex_fine;i++)
    {
      for(j=0;j<3;j++)
	{
	  if(fine_obj->myvertexs[i].location(j)>loc_max(j))
	    {
	      loc_max(j)=fine_obj->myvertexs[i].location(j);
	    }
	  if(fine_obj->myvertexs[i].location(j)<loc_min(j))
	    {
	      loc_min(j)=fine_obj->myvertexs[i].location(j);
	    }
	}
    }

  myvector diag=loc_max-loc_min; double norm=diag.len();

  myvector max_error=fine_obj->myvertexs[0].location-coarsen_obj->myvertexs[0].location;
  myvector error_now;  double error_all=0;
  for(i=0;i<num_vertex_coarsen;i++)
    {
      error_now=fine_obj->myvertexs[i].location-coarsen_obj->myvertexs[i].location;
      error_all+=error_now.len();
      if(error_now.len()>max_error.len())
	{
	  max_error=error_now;
	}
    }

  double average_error_value=error_all/(double)num_vertex_coarsen/norm;
  double max_error_value=max_error.len()/norm;

  myio.save_ERROR_TXT(average_error_value,average_error_txt);
  myio.save_ERROR_TXT(max_error_value,max_error_txt);
  
  delete fine_obj;
  delete coarsen_obj;  
}
