#ifndef _MATE_MAKER_
#define _MATE_MAKER_

#include "head.h"
#include "hexahedron.h"
class mate_maker
{
 public:
  ~mate_maker(){}
  mate_maker(){}
  mate_maker(boost::property_tree::ptree &para_tree);
  int cal_mate(hexahedron &hex_now,double y,double p);
};

#endif
