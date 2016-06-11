 #ifndef __DRIVER_H
#define __DRIVER_H

#include <vector>
#include <list>
#include <utility> 

using namespace std;

class Driver {
 private:
  unsigned int index;
  unsigned int maxDrivingDuration;
  vector<pair <int,int> > timeWindows;
  unsigned int trailer;
  unsigned int minInterSHIFTDURATION;
  double timeCost;
  
public:
  Driver();
  unsigned int getIndex();
  void setIndex(unsigned int index);
  unsigned int getMaxDrivingDuration();
  void setMaxDrivingDuration(unsigned int maxDD);
  vector<pair <int,int> > getTimeWindows();
  void setTimeWindows(vector<pair <int,int> > tw);
  unsigned int getTrailer();
  void setTrailer(unsigned int t);
  unsigned int getMinInterSHIFTDURATION();
  void setMinInterSHIFTDURATION(unsigned int t);
  double getTimeCost();
  void setTimeCost(double t);
  

};

#endif
