#include "face.h"

face::face(const size_t (&index_vertex)[4])
{
  size_t i;
  for(i=0;i<4;++i)
    {
      this->index_vertex[i]=index_vertex[i];
    }
  return;
}
