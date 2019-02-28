#include "exchanger.h"
#include "io.h"
#include "object.h"
using namespace std;
exchanger::exchanger(boost::property_tree::ptree &para_tree)
{
  string file_vtk_fine=para_tree.get<string >("file_vtk_fine.value");
  string mat_in_fine=para_tree.get<string >("mat_in_fine.value");
  string mat_out_fine=para_tree.get<string>("mat_out_fine.value");
  
  io myio=io();
  object* fine_obj=new object();
  myio.getVertexAndHex(fine_obj->myvertexs,fine_obj->myhexahedrons,file_vtk_fine);
  fine_obj->num_hexahedrons=fine_obj->myhexahedrons.size();
  fine_obj->num_vertex=fine_obj->myvertexs.size();
  myio.getMatParaFromTXT(fine_obj,mat_in_fine);
  myio.saveMatPara(fine_obj,mat_out_fine,1);
  delete fine_obj;
  
}
