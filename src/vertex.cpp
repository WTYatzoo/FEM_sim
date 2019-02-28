#include "vertex.h"
using namespace std;

vertex::vertex()
{
  ;
}

vertex::vertex(const myvector &location)
{
  this->location=location;
  this->location_original=location;
  this->location_maybe=location;
  this->velocity=myvector(0,0,0);
  this->force_external=myvector(0,0,0);
  isFixed=0; // all free when at beginning
}

vertex::~vertex()
{
  ;
}
