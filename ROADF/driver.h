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

/*
                 if(f==0)
                    solution.tankQuantities[c][f] -= this->customers[c].getForecast()[f];
                else
                    solution.tankQuantities[c][f] += solution.tankQuantities[c][f-1] - this->customers[c].getForecast()[f];
                    */
/*
    for(int s=0; s<solution.getShifts().size(); s++){
        for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++){
            Operation op = solution.getShifts()[s].getOperations()[o];
            if(op.getPoint() > 0){
                for(int f=op.getArrival()+this->customers[op.getPoint()].getSetupTime(); f<this->customers[op.getPoint()].getForecast().size()*60; f++){
                    solution.tankQuantities[op.getPoint()][f] += op.getQuantity();
                    solution.trailerQuantities[op.getPoint()][f] -= op.getQuantity();
                }
            }
        }
    }
    cout<<"\n";
    */
#endif
