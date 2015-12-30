#ifndef __LOCATION_H
#define __LOCATION_H

#include<vector>

#include "trailer.h"


using namespace std;

class Customer{
private:
  unsigned int index;
  unsigned int setupTime;
  vector<unsigned int> allowedTrailers;
  vector<double> forecast;
  double capacity;
  double initialTankQuantity;
  double safetyLevel;
    
public:
  Customer();
  unsigned int getIndex();
  void setIndex(unsigned int ind);
  unsigned int getSetupTime();
  void setSetupTime(unsigned int st);
  vector<unsigned int> getAllowedTrailers();
  void setAllowedTrailers(vector<unsigned int> at);
  vector<double> getForecast();
  void setForecast(vector<double> f);
  double getCapacity();
  void setCapacity(double c);
  double getInitialTankQuantity();
  void setInitialTankQuantity(double itq);
  double getSafetyLevel();
  void setSafetyLevel(double sl);
};
#endif
