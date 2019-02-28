#ifndef _FACE_
#define _FACE_
#include "myvector.h"

class face
{
 public:
  size_t index_vertex[4];  //这个面的顶点索引是从小到大排列的 
  myvector normal_ori; //指向外面
  size_t num_hex; //拥有这个面的四面体的数目

  face(){}
  ~face(){}
  face(const size_t (&index_vertex)[4]);
};
#endif
