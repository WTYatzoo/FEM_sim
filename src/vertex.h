#ifndef _VERTEX_
#define _VERTEX_

#include "head.h"
#include "myvector.h"

class vertex
{
 public:
  myvector location;
  myvector velocity;
  myvector location_original;
  myvector location_lastFrame;
  myvector velocity_lastFrame;
  myvector location_maybe;
  myvector force_external;
  double mass;
  double scalar_field_value; // for gyroid structure based material generation 
  int isFixed; //  1:fixed points¡¡0:free points
  vertex();
  vertex(const myvector &location);
  ~vertex();
};

#endif
