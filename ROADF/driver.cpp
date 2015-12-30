#ifndef __DRIVER_CPP
#define __DRIVER_CPP

#include "driver.h"

  Driver::Driver(){}
  unsigned int Driver::getIndex()    {return this->index;}
  void Driver::setIndex(unsigned int ind)	{this->index=ind;}
  unsigned int Driver::getMaxDrivingDuration()    {return this->maxDrivingDuration;}
  void Driver::setMaxDrivingDuration(unsigned int maxDD)	{this->maxDrivingDuration=maxDD;}
  std::vector<std::pair <int,int> > Driver::getTimeWindows()	{return this->timeWindows;}
  void Driver::setTimeWindows(std::vector<std::pair <int,int> > tw)	{this->timeWindows = tw;}
  unsigned int Driver::getTrailer()    {return this->trailer;}
  void Driver::setTrailer(unsigned int t)	{this->trailer=t;}
  unsigned int Driver::getMinInterSHIFTDURATION()    {return this->minInterSHIFTDURATION;}
  void Driver::setMinInterSHIFTDURATION(unsigned int s)	{this->minInterSHIFTDURATION=s;}
  double Driver::getTimeCost()    {return this->timeCost;}
  void Driver::setTimeCost(double t)	{this->timeCost=t;}
  
  
#endif
