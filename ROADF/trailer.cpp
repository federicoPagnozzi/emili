#ifndef __TRAILER_CPP
#define __TRAILER_CPP

#include "trailer.h"

  Trailer::Trailer(){}
  unsigned int Trailer::getIndex()    {return this->index;}
  void Trailer::setIndex(unsigned int ind)	{this->index=ind;}
  double Trailer::getCapacity()    {return this->capacity;}
  void Trailer::setCapacity(double c)	{this->capacity=c;}
  double Trailer::getInitialQuantity()    {return this->initialQuantity;}
  void Trailer::setInitialQuantity(double iq)	{this->initialQuantity=iq;}
  double Trailer::getDistanceCost()    {return this->distanceCost;}
  void Trailer::setDistanceCost(double dc)	{this->distanceCost=dc;}
  
  
#endif
