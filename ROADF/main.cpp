#include <iostream>
#include <sstream>

#include "tinyxml.h"
#include "tinystr.h"
#include "matrix.h"
#include "instance.h"
#include "solution.h"
#include "driver.h"

#include <string>
#include <vector>

int maino()
{
  srand (time(NULL));
   
  Instance instance;
  irpSolution solution;
  instance.loadInstance("Instance_V_1.5.xml");
  solution.loadSolution("Solution_V_1.5.xml");

 /* if(instance.checkFeasibility(solution)){
    cout<<"NOT FEASIBLE!";
    exit(0);
  }*/
  
  
  double par = 1.0;
  /*
  Solution rs=instance.randomSolution(par);
  
  if(instance.checkFeasibility(rs)){
//    rs=instance.randomSolution(par);
    cout<<"NOT FEASIBLE! "<<par<<"\n";
//    par += 0.01;
    rs.saveSolution("MySolution_V_1.1.xml");
  }
  else
    rs.saveSolution("MySolution_V_1.3.xml");*/
  


//  solution = instance.backTrackingRandomSolution(1.0, 1.0, 1);
    cout<<"LR: "<<instance.computeObjective(solution)<<"\n";
  if(instance.checkFeasibility(solution, false))
    cout<<"\nNOT FEASIBLE! "<<par<<"\n";
  solution.saveSolution("BackSolution_V_1.1.xml");
}
