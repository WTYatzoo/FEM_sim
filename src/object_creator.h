#ifndef _OBJECT_CREATOR_
#define _OBJECT_CREATOR_

#include "head.h"
#include "myvector.h"
#include "hexahedron.h"
#include "object.h"

class object_creator //create a pair of coarsen & fine meshes for testing numerical coarsening, which need 1.粗网格的顶点序号和细网格对应顶点序号相同 2.细网格每8个对应粗网格一个
{
 public:
  double length_fine,width_fine,height_fine; 
  int*** index_for_vertex; //因为需要赋值为-1,所以用int 
  size_t lc_fine,wc_fine,hc_fine; //节点数
  int level; //for the edge: 1: 2 or 1: 4 or 1:8
  double dmetric_fine;
  double dmetric_coarsen;
  double PoissonRatio,YoungModulus1,YoungModulus2;
  std::string object_name;
  std::string out_dir;
  object_creator();
  ~object_creator();
  object_creator(boost::property_tree::ptree &para_tree);
  int create_object();
  int create_object_bridge();
  int create_object_bridge_sin();
  int create_fine_object();
  int bind_mat(hexahedron &hexahedron_here,const double &PoissonRatio,const double &YoungModulus);
};

#endif
