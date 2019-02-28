#include "head.h"
#include "object_creator.h"
#include "simulator.h"
#include "harmonic_solver.h"
#include "coarsener.h"
#include "exchanger.h"
#include "spectrum_analyser.h"
#include "CtoFer.h"
#include "mate_maker.h"
#include "mate_factory.h"
#include "error_measurer.h"
using namespace std;

void readcmdline(int argc, char* argv[],boost::property_tree::ptree &para_tree)
{
  size_t i;
  for(i=1;i<argc;++i)
    {
      string para_here=argv[i];
      size_t pos=para_here.find("=");
      if(pos!= string::npos)
	{
	  string key=para_here.substr(0,pos);
	  string value=para_here.substr(pos+1);
	  para_tree.put(key+".value",value);
	  printf("--[cmdline para] %s %s \n",key.c_str(),value.c_str());
	}
    }
  return;
}
int main(int argc, char *argv[])
{
  // test access point
  boost::property_tree::ptree para_tree;
  readcmdline(argc,argv,para_tree);
  string prog=para_tree.get<string>("prog.value");
  if(prog=="object_creator")
    {
      object_creator* myobject_creator=new object_creator(para_tree);
      delete myobject_creator;
    }
  else if(prog=="harmonic_solver")
    {
      harmonic_solver* myharmonic_slover=new harmonic_solver(para_tree);
      delete myharmonic_slover;
    }
  else if(prog=="simulator")
    {
      simulator* mysimulator=new simulator(para_tree);
      delete mysimulator;
    }
  else if(prog=="coarsener")
    {
      coarsener* mycoarsener=new coarsener(para_tree);
      delete mycoarsener;
    }
  else if(prog=="exchanger")
    {
      exchanger* myexchanger=new exchanger(para_tree);
      delete myexchanger;
    }
  else if(prog=="spectrum_analyser")
    {
      spectrum_analyser* myspectrum_analyser=new spectrum_analyser(para_tree);
      delete myspectrum_analyser;
    }
  else if(prog=="CtoFer")
    {
      CtoFer* myCtoFer=new CtoFer(para_tree);
      delete myCtoFer;
    }
  else if(prog=="mate_maker")
    {
      mate_maker* mymate_maker=new mate_maker(para_tree);
      delete mymate_maker;
    }
  else if(prog=="mate_factory")
    {
      mate_factory* mymate_factory=new mate_factory(para_tree);
      delete mymate_factory;
    }
  else if(prog=="error_measurer")
    {
      error_measurer* myerror_measurer=new error_measurer(para_tree);
      delete myerror_measurer;
    }
  return 0;
}
