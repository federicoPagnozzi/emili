#ifndef __TRAILER_H
#define __TRAILER_H

#include <vector>
#include <list>
#include <utility> 

class Trailer {
 private:
  unsigned int index;
  double capacity;
  double initialQuantity;
  double distanceCost;
  
public:
  Trailer();
  unsigned int getIndex();
  void setIndex(unsigned int index);
  double getCapacity();
  void setCapacity(double c);
  double getInitialQuantity();
  void setInitialQuantity(double ic);
  double getDistanceCost();
  void setDistanceCost(double dc);

};

#endif
