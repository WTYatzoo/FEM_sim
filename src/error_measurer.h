#ifndef _ERROR_MEASURER_
#define _ERROR_MEASURER_

#include "head.h"

class error_measurer
{
 public:
  ~error_measurer(){}
  error_measurer(){}
  error_measurer(boost::property_tree::ptree &para_tree);
};

#endif
