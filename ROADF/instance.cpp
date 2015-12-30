#ifndef __INSTANCE_CPP
#define __INSTANCE_CPP

#include "instance.h"
#include "tinyxml.h"
#include "tinystr.h"

#include <cmath>
#include <iomanip>

#define EPSILON 0.000001

Instance::Instance()
{
}

int Instance::getUnit()	{return this->unit;}
int Instance::getHorizon()	{return this->horizon;}
vector<vector<double> > Instance::getTimeMatrices()	{return this->timeMatrices;}
vector<Driver> Instance::getDrivers()	{return this->drivers;}
vector<Trailer> Instance::getTrailers()	{return this->trailers;}
vector<Customer> Instance::getCustomers()	{return this->customers;}
vector<vector<double> > Instance::getDistMatrices()	{return this->distMatrices;}


/**
 * Read xml instance file
 */
void Instance::loadInstance(const char* pFilename)
{
      int index;
      double d;
      TiXmlDocument doc(pFilename);

      if (!doc.LoadFile()) return;
      
      TiXmlHandle hDoc(&doc);
      TiXmlElement* pInnerElem;
      TiXmlHandle hRoot(0);


      pInnerElem=hDoc.FirstChildElement().Element();
      // should always have a valid root but handle gracefully if it does
      if (!pInnerElem) return;
      name=pInnerElem->Value();

      // save this for later
      hRoot=TiXmlHandle(pInnerElem);
      
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild().Element();
      
      //Read unit
      istringstream ( pInnerElem->GetText()) >> this->unit;
//    cout<<pInnerElem->GetText()<<" "<<this->unit<<"\n";
      
      //Read horizon
      pInnerElem=pInnerElem->NextSiblingElement();
      istringstream ( pInnerElem->GetText()) >> this->horizon;
//    cout<<pInnerElem->Value()<<" "<<this->horizon<<"\n";
      
      //Read time matrix
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("timeMatrices").FirstChild().Element();
      int i = 0;
      for( pInnerElem; pInnerElem; pInnerElem=pInnerElem->NextSiblingElement())
      {
	      int j = 0;
	      vector<double> vec;
	      this->timeMatrices.push_back(vec);
	      TiXmlElement *pMatrixElem = pInnerElem->FirstChild()->ToElement();
	      for( pMatrixElem ; pMatrixElem;pMatrixElem =pMatrixElem ->NextSiblingElement())
	      {
		int elem;
		istringstream(pMatrixElem ->GetText()) >> elem;
		this->timeMatrices.at(i).push_back(elem);
//		cout<<this->timeMatrices.at(i).at(j)<<" ";
		j++;
	      }
//	      cout<<"\n";
	      i++;
      }
      
      //Read drivers data
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("drivers").FirstChild().Element();
      while(pInnerElem!= NULL and strcmp(pInnerElem->Value(), "IRP_Roadef_Challenge_Instance_driver")==0)
      {
	Driver d;
	
	TiXmlElement *pDriverElem = pInnerElem->FirstChild()->ToElement();
	
	//Read driver's index
	istringstream(pDriverElem->GetText()) >> index;
	d.setIndex(index);
//	cout<<pDriverElem->Value()<< " " << d.getIndex()<<"\n";
	
	//Read driver's masimum driving duration
	pDriverElem=pDriverElem->NextSiblingElement();	
	istringstream(pDriverElem->GetText()) >> index;
	d.setMaxDrivingDuration(index);	
//	cout<<pDriverElem->Value()<< " " << d.getMaxDrivingDuration()<<"\n";
	
	//Read driver's time windows (start,end)
	pDriverElem=pDriverElem->NextSiblingElement();
	TiXmlElement *pTimeWindowElem=pDriverElem->FirstChild()->ToElement();
	vector<pair<int,int> > tw;
	while(pTimeWindowElem!= NULL and strcmp(pTimeWindowElem->Value(), "TimeWindow")==0){
	   pair <int,int> p (0,0);
	   istringstream(pTimeWindowElem->FirstChild()->ToElement()->GetText()) >> p.first;
	   istringstream(pTimeWindowElem->FirstChild()->ToElement()->NextSiblingElement()->GetText()) >> p.second;
	   tw.push_back(p);
	   pTimeWindowElem=pTimeWindowElem->NextSiblingElement();
	}
	d.setTimeWindows(tw);
	
	//Read driver's trailer
	pDriverElem=pDriverElem->NextSiblingElement();	
	istringstream(pDriverElem->GetText()) >> index;
	d.setTrailer(index);	
//	cout<<pDriverElem->Value()<< " " << d.getTrailer()<<"\n";
	
	//Read driver's minimum interval shift duration
	pDriverElem=pDriverElem->NextSiblingElement();	
	istringstream(pDriverElem->GetText()) >> index;
	d.setMinInterSHIFTDURATION(index);	
//	cout<<pDriverElem->Value()<< " " << d.getMinInterSHIFTDURATION()<<"\n";
	
	//Read driver's time cost
	double tc;
	pDriverElem=pDriverElem->NextSiblingElement();	
	istringstream(pDriverElem->GetText()) >> tc;
	d.setTimeCost(tc);	
//	cout<<pDriverElem->Value()<< " " <<d.getTimeCost()<<"\n";
	
	this->drivers.push_back(d);
	
	pInnerElem=pInnerElem->NextSiblingElement();
      }
      
      //Read trailers data 
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("trailers").FirstChild().Element();
      while(pInnerElem!= NULL and strcmp(pInnerElem->Value(), "IRP_Roadef_Challenge_Instance_Trailers")==0)
      {
	Trailer t;
	TiXmlElement *pTrailerElem = pInnerElem->FirstChild()->ToElement();
	
	//Read trailer's index	
	istringstream(pTrailerElem->GetText()) >> index;
	t.setIndex(index);	
//	cout<<pTrailerElem->Value()<< " " << t.getIndex()<<"\n";
	
	//Read trailer's capacity
	pTrailerElem=pTrailerElem->NextSiblingElement();	
	istringstream(pTrailerElem->GetText()) >> d;
	t.setCapacity(d);	
//	cout<<pTrailerElem->Value()<< " " << t.getCapacity()<<"\n";
	
	//Read trailer's initial quantity
	pTrailerElem=pTrailerElem->NextSiblingElement();	
	istringstream(pTrailerElem->GetText()) >> d;
	t.setInitialQuantity(d);	
//	cout<<pTrailerElem->Value()<< " " << t.getInitialQuantity()<<"\n";
	
	//Read trailer's distance cost
	pTrailerElem=pTrailerElem->NextSiblingElement();	
	istringstream(pTrailerElem->GetText()) >> d;
	t.setDistanceCost(d);	
//	cout<<pTrailerElem->Value()<< " " << t.getDistanceCost()<<"\n";
	
	this->trailers.push_back(t);
	
	pInnerElem=pInnerElem->NextSiblingElement();
      }
      
      //Read bases data
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("bases").FirstChild().Element();
      while(pInnerElem!= NULL and strcmp(pInnerElem->Value(), "index")==0)
      {
	Customer b;
	int index;

	//Read base's index
	istringstream(pInnerElem->GetText()) >> index;
	b.setIndex(index);	
//	cout<<pInnerElem->Value()<< " " << b.getIndex()<<"\n";

	/*this->bases.push_back(b);*/this->customers.push_back(b);
	
	pInnerElem=pInnerElem->NextSiblingElement();
      }
      
      //Read sources data
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("sources").FirstChild().Element();
      while(pInnerElem!= NULL and strcmp(pInnerElem->Value(), "IRP_Roadef_Challenge_Instance_Sources")==0)
      {
	Customer s;
	double d;
	TiXmlElement *pSourceElem = pInnerElem->FirstChild()->ToElement();
	
	//Read source's index
	istringstream(pSourceElem->GetText()) >> index;
	s.setIndex(index);	
//	cout<<pSourceElem->Value()<< " " << s.getIndex()<<"\n";
	
	//Read source's setup time
	pSourceElem=pSourceElem->NextSiblingElement();	
	istringstream(pSourceElem->GetText()) >> index;
	s.setSetupTime(index);	
//	cout<<pSourceElem->Value()<< " " << s.getSetupTime()<<"\n";
	
	/*this->sources.push_back(s);*/this->customers.push_back(s);
	
	pInnerElem=pInnerElem->NextSiblingElement();
      }
      
      //Read customers data
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("customers").FirstChild().Element();
      while(pInnerElem!= NULL and strcmp(pInnerElem->Value(), "IRP_Roadef_Challenge_Instance_Customers")==0)
      {
	Customer c;
	
	TiXmlElement *pCustomerElem = pInnerElem->FirstChild()->ToElement();
	
	//Read customer's indes
	istringstream(pCustomerElem->GetText()) >> index;
	c.setIndex(index);
//	cout<<pCustomerElem->Value()<< " " << c.getIndex()<<"\n";
	
	//Read customer's setup time
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> index;
	c.setSetupTime(index);	
//	cout<<pCustomerElem->Value()<< " " << c.getSetupTime()<<"\n";
	
        //Read customer's allowed trailers
	pCustomerElem=pCustomerElem->NextSiblingElement();
	TiXmlElement *pAllowedTrailersElem=pCustomerElem->FirstChild()->ToElement();
	vector<unsigned int> at;
	while(pAllowedTrailersElem!= NULL and strcmp(pAllowedTrailersElem->Value(), "int")==0){
	    
	    istringstream(pAllowedTrailersElem->ToElement()->GetText()) >> index;
	    at.push_back(index);
//	    cout<<pAllowedTrailersElem->Value()<< " " << at.at(0)<<"\n";
	    pAllowedTrailersElem=pAllowedTrailersElem->NextSiblingElement();
	}
	c.setAllowedTrailers(at);
	
	//Read customer's forecast
	pCustomerElem=pCustomerElem->NextSiblingElement();
	TiXmlElement *pForecastElem=pCustomerElem->FirstChild()->ToElement();
	vector<double> f;
	while(pForecastElem!= NULL and strcmp(pForecastElem->Value(), "double")==0){
	    istringstream(pForecastElem->ToElement()->GetText()) >> d;
	    f.push_back(d);
//	    cout<<"here here"<<pForecastElem->Value()<< " " << f.at(0)<<"\n";
	    pForecastElem=pForecastElem->NextSiblingElement();
	}
	c.setForecast(f);
	
	//Read customer's capacity
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> d;
	c.setCapacity(d);	
//	cout<<pCustomerElem->Value()<< " " << c.getCapacity()<<"\n";
	
	//Read customer's initial tank quantity
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> d;
	c.setInitialTankQuantity(d);	
//	cout<<pCustomerElem->Value()<< " " << c.getInitialTankQuantity()<<"\n";
	
	//Read customer's safety level
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> d;
	c.setSafetyLevel(d);	
//	cout<<pCustomerElem->Value()<< " " << c.getSafetyLevel()<<"\n";
	 
	this->customers.push_back(c);
	
	pInnerElem=pInnerElem->NextSiblingElement();
      }
      
      //Read distance matrix
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("DistMatrices").FirstChild().Element();
      i=0;
      for( pInnerElem; pInnerElem; pInnerElem=pInnerElem->NextSiblingElement())
      {
	      int j = 0;
	      vector<double> vec;
	      this->distMatrices.push_back(vec);
//	      cout<<pInnerElem ->Value()<<": ";
	      TiXmlElement *pMatrixElem = pInnerElem->FirstChild()->ToElement();
	      for( pMatrixElem ; pMatrixElem;pMatrixElem =pMatrixElem ->NextSiblingElement())
	      {
		int elem;
		istringstream(pMatrixElem ->GetText()) >> elem;
		this->distMatrices.at(i).push_back(elem);
	//	cout<<this->distMatrices.at(i).at(j)<<" ";
		j++;
	      }
	//      cout<<"\n";
	      i++;
      }

}



/**
 * Driver Constraint:
 * [DRI01 | Inter-Shifts duration]
 */
/*
bool Instance::dri01(irpSolution solution){
  unsigned int end1,end2;
  for(int d=0; d<this->drivers.size(); d++){
      
    Driver driver=this->drivers[d];
    for(int s1=0; s1<solution.getShifts().size()-1; s1++){
      for(int s2=s1+1; s2<solution.getShifts().size(); s2++){
	
	Shift shift1=solution.getShifts()[s1];
	Shift shift2=solution.getShifts()[s2];
    if(shift1.getDriver()==driver.getIndex()  && shift2.getDriver()==driver.getIndex()){
	  
	  end1 = shift1.getOperations().back().getArrival() + this->customers.at(shift1.getOperations().back().getPoint()).getSetupTime()
	         + this->timeMatrices[shift1.getOperations().back().getPoint()][0];
	  end2 = shift2.getOperations().back().getArrival() + this->customers.at(shift2.getOperations().back().getPoint()).getSetupTime()
		 + this->timeMatrices[shift2.getOperations().back().getPoint()][0];
//	  cout<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<shift1.getStart()<<" "<<end1<<" \n";
	  if(  not(  shift2.getStart() >= end1  + driver.getMinInterSHIFTDURATION() 
          //  or shift1.getStart()  >= end2 + driver.getMinInterSHIFTDURATION()
                 )
           ){
	    cout<<"dri01 NOT FEASIBLE!!!!\n";
//	    cout<<shift1.getStart()<<" "<<end1<<" "<<shift2.getStart()<<" "<<end2<<"\n";
        //exit(0);
	    return true;
	  }
	}
	
      }
    }
  }
    return false;
}
*/

bool Instance::dri01(irpSolution solution){
  unsigned int end1,end2;

    for(int s1=0; s1<solution.getShifts().size()-1; s1++){
      for(int s2=s1+1; s2<solution.getShifts().size(); s2++){

        Shift shift1 = solution.getShifts()[s1];
        Shift shift2 = solution.getShifts()[s2];

        if(shift1.getDriver() == shift2.getDriver()){
          Driver driver = this->drivers[shift1.getDriver()];
          end1 = shift1.getOperations().back().getArrival() + this->customers.at(shift1.getOperations().back().getPoint()).getSetupTime()
                 + this->timeMatrices[shift1.getOperations().back().getPoint()][0];
          end2 = shift2.getOperations().back().getArrival() + this->customers.at(shift2.getOperations().back().getPoint()).getSetupTime()
             + this->timeMatrices[shift2.getOperations().back().getPoint()][0];
    //	  cout<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<shift1.getStart()<<" "<<end1<<" \n";
          if(  not(  shift2.getStart() >= end1  + driver.getMinInterSHIFTDURATION()
                /*or shift1.getStart()  >= end2 + driver.getMinInterSHIFTDURATION() */)
               ){
            cout<<"dri01 NOT FEASIBLE!!!!\n";
    //	    cout<<shift1.getStart()<<" "<<end1<<" "<<shift2.getStart()<<" "<<end2<<"\n";
            //exit(0);
            return true;
          }
        }

          }
    }
    return false;
}

/**
 * Driver Constraint:
 * [DRI03 | Respect of maximal driving time]
 */
bool Instance::dri03(irpSolution solution){
  
  for(int s=0; s<solution.getShifts().size(); s++){
    unsigned int cumulatedDrivingTime = 0;
    
// cout<<"CHECK SHIFT: "<<  solution.getShifts().at(s).getIndex()<<" ";
   vector<Operation> operations = solution.getShifts()[s].getOperations();
   cumulatedDrivingTime += this->timeMatrices[0][operations.front().getPoint()];
   for(int o=1; o<operations.size(); o++){
     
     cumulatedDrivingTime+=this->timeMatrices[operations[o-1].getPoint()][operations[o].getPoint()];
//     cout<<operations.at(o-1).getPoint()<<" "<<operations.at(o).getPoint()<<" "<<cumulatedDrivingTime<<"\n";  
   }
   cumulatedDrivingTime += this->timeMatrices[operations.back().getPoint()][0];
   if(not(cumulatedDrivingTime <= this->drivers.at(solution.getShifts()[s].getDriver()).getMaxDrivingDuration())){
     cout<<"DRI03 NOT FEASIBLE!!!!";
     /* exit(0);*/
      return true;
   }
//    cout<<"\n";
  }
  return false; 
}

/**
 * Driver Constraint:
 * [DRI08 | Time windows of the drivers]
 */
bool Instance::dri08(irpSolution solution){
  
  for(int s=0; s<solution.getShifts().size(); s++){
    
    bool flag=true;
    Shift shift=solution.getShifts()[s];
    vector<pair<int,int> > timeWindows=this->drivers[shift.getDriver()].getTimeWindows();
//    cout<<"SHIFT: "<<shift.getIndex()<<"\n";
    for(int tw=0; tw<timeWindows.size(); tw++){
      
      unsigned int end = shift.getOperations().back().getArrival() + this->customers.at(shift.getOperations().back().getPoint()).getSetupTime()
			 + this->timeMatrices[shift.getOperations().back().getPoint()][0];
//      cout<<timeWindows[tw].first<<" "<<timeWindows[tw].second<<" "<<shift.getStart()<<" "<<end<<"\n";
      if(shift.getStart() >= timeWindows[tw].first && end <= timeWindows[tw].second
              /*&& shift.getStart() <= timeWindows[tw].second && end >= timeWindows[tw].first*/
        && shift.getStart()<= end && timeWindows[tw].first<=timeWindows[tw].second){
            flag = false;
      }

    }
    
    if(flag){
       cout<<"DRI08 NOT FEASIBLE!!!!";
       
	/*exit(0);*/
	return true;
    }
  }
  return false; 
}

	
/**
 * Trailer Constraint:
 * [TL01 | Different shifts of the same trailer cannot overlap in time]
 */
/*
bool Instance::tl01(irpSolution solution){
  unsigned int end1,end2;
  for(int t=0; t<this->trailers.size(); t++){
    
    Trailer trailer = this->trailers[t];
    for(int s1=0; s1<solution.getShifts().size(); s1++){
      for(int s2=s1+1; s2<solution.getShifts().size(); s2++){

    Shift shift1 = solution.getShifts()[s1];
    Shift shift2 = solution.getShifts()[s2];
    if(shift1.getTrailer()==trailer.getIndex() && shift2.getTrailer()==trailer.getIndex()){
	  end1 = shift1.getOperations().back().getArrival() + this->customers.at(shift1.getOperations().back().getPoint()).getSetupTime()
		  + this->timeMatrices[shift1.getOperations().back().getPoint()][0];
	  end2 = shift2.getOperations().back().getArrival() + this->customers.at(shift2.getOperations().back().getPoint()).getSetupTime()
		  + this->timeMatrices[shift2.getOperations().back().getPoint()][0];
//	  cout<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<end1<<" \n";

      if(  not(  shift2.getStart() >= end1 //or shift1.getStart()  >= end2
                 )){
	     cout<<"TL01 NOT FEASIBLE!!!!";
       //exit(0);
	    return true;
	  }

	}
	
      }
    }
    
  }  
  return false; 
}
*/
bool Instance::tl01(irpSolution solution){
  unsigned int end1,end2;

    for(int s1=0; s1<solution.getShifts().size(); s1++){
      for(int s2=s1+1; s2<solution.getShifts().size(); s2++){

        Shift shift1 = solution.getShifts()[s1];
        Shift shift2 = solution.getShifts()[s2];
        if(shift1.getTrailer() == shift2.getTrailer()){
          end1 = shift1.getOperations().back().getArrival() + this->customers.at(shift1.getOperations().back().getPoint()).getSetupTime()
              + this->timeMatrices[shift1.getOperations().back().getPoint()][0];
          end2 = shift2.getOperations().back().getArrival() + this->customers.at(shift2.getOperations().back().getPoint()).getSetupTime()
              + this->timeMatrices[shift2.getOperations().back().getPoint()][0];
    //	  cout<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<end1<<" \n";

          if(  not(  shift2.getStart() >= end1 //or shift1.getStart()  >= end2
                     )){
             cout<<"TL01 NOT FEASIBLE!!!!";
           //exit(0);
            return true;
          }

        }

      }
    }

  return false;
}


/**
 * Trailer Constraint:
 * [TL03 | The trailer attached to a driver in a shift must be compatible]
 */
bool Instance::tl03(irpSolution solution){
  
  for(int s=0; s<solution.getShifts().size(); s++){
    
//    cout<<solution.getShifts()[s].getTrailer().getIndex()<<" "<<solution.getShifts()[s].getDriver().getIndex()<<" "
//    <<this->drivers[solution.getShifts()[s].getDriver().getIndex()].getTrailer()<<"\n";
    if( not(solution.getShifts()[s].getTrailer()==this->drivers[solution.getShifts()[s].getDriver()].getTrailer() ) ){
       cout<<"TL03 NOT FEASIBLE!!!!";
	/*exit(0);*/
	return true;
    }
  }
  return false; 
}

/**
 * Site Constraint:
 * [DYN01 | Respect of tank capacity for each site]
 * [SHI11 | Some product must be loaded or delivered]
 * [QS02  | Run-out avoidance]
 */
double Instance::dyn01(irpSolution solution, bool feasibility){
  
  double unfeasibilityCounter = 0.0;

  vector<vector<double> > tankQuantities(this->customers.size());
  for(int p=0; p<this->customers.size();p++){
    if(p>1){
      vector<double> horizon(this->customers[p].getForecast().size()*60,0.0);
      tankQuantities[p]=horizon;
      tankQuantities[p][0] = this->customers[p].getInitialTankQuantity();
    }
  }
  
  for(int s=0; s<solution.getShifts().size(); s++){
    Shift shift=solution.getShifts()[s];
    vector<Operation> operations = shift.getOperations();
    for(int o=0; o<operations.size(); o++){

        if(operations[o].getPoint() > 1){
          if(not (operations[o].getQuantity() >= -EPSILON)){
//             cout<<"SHI11 NOT FEASIBLE!!!!";
                cout<<"\n"<<operations[o].getPoint()<<" "<<operations[o].getQuantity()<<" ";
                unfeasibilityCounter = operations[o].getArrival();
              /*exit(0);*/
          //    return true;
          }
          tankQuantities[operations[o].getPoint()][(int)operations[o].getArrival()] += operations[o].getQuantity();
        }
        else if(operations[o].getPoint() == 1)
          if(not (operations[o].getQuantity() <= EPSILON)){
//             cout<<"SHI11 NOT FEASIBLE!!!!";
             unfeasibilityCounter = operations[o].getArrival();
              /*exit(0);*/
          //    return true;
          }
    }
  }  

  double customerCounter = 1;
  bool customerFlag = false;
  for(int p=0; p<this->customers.size(); p++){
    if(p>1){
      for(int f=0; f<this->customers[p].getForecast().size()*60; f++){
        if(f%60 == 0)
        tankQuantities[p][f] -= this->customers[p].getForecast()[(int)f/60];
        if(f>=1)
          tankQuantities[p][f] += tankQuantities[p][f-1];	  

        if( not(  (this->customers[p].getCapacity() - tankQuantities[p][f] >= -EPSILON or abs(this->customers[p].getCapacity() - tankQuantities[p][f]) <= EPSILON)
          && tankQuantities[p][f] >= -EPSILON
          && ( tankQuantities[p][f] >= this->customers[p].getSafetyLevel() or abs(tankQuantities[p][f] - this->customers[p].getSafetyLevel()) <= EPSILON)
                   )
             ){
//               cout<<"DYN01 QS02 NOT FEASIBLE!!!!";
//               cout<<"\n"<<p<<" "<<f<<" "<<tankQuantities[p][f]<<" "<<
//            this->customers[p].getCapacity()<<" "<<this->customers[p].getSafetyLevel()<<"\n";
               if(not(customerFlag)){
                  unfeasibilityCounter += abs(this->horizon*60 - f);
                  customerFlag = true;
               }
            //   exit(0);
            //    return true;
            }

      }
    }
//    if(customerFlag)
//        customerCounter++;
    customerFlag = false;
  }
//  return false;
  return unfeasibilityCounter;
}

/**
 * Shift Constraint:
 * [SHI02 | Arrival at point requires travelling time from previous point]
 * [SHI03 | Loading and delivery operations take a constant time]
 */
bool Instance::shi02(irpSolution solution){
  
    for(int s=0; s<solution.getShifts().size(); s++){
        Shift shift=solution.getShifts()[s];
        vector<Operation> operations = shift.getOperations();

        double departure=shift.getStart();
        if(not(  operations[0].getArrival() >= departure + this->timeMatrices[0][operations[0].getPoint()] )){
        cout<<"SHI03 NOT FEASIBLE!!!!";
         /*exit(0);*/
         return true;
        }

        for(int o=1; o<operations.size(); o++){
          departure = operations[o-1].getArrival() + this->customers[operations[o-1].getPoint()].getSetupTime();
    //      if(s==1)
    //      cout<<operations[o].getArrival()<<" "<<departure + this->timeMatrices[operations[o-1].getPoint()][operations[o].getPoint()]<<"\n";
          if(not(  operations[o].getArrival() >= departure + this->timeMatrices[operations[o-1].getPoint()][operations[o].getPoint()] )){
         cout<<"SHI02 NOT FEASIBLE!!!!";
         /*exit(0);*/
         return true;
          }
        }
  } 
  return false; 
}

/**
 * Shift Constraint:
 * [SHI05 | Delivery operations require the customer site to be accessible for trailer]
 */
bool Instance::shi05(irpSolution solution){
  
   for(int s=0; s<solution.getShifts().size(); s++){
    Shift shift=solution.getShifts()[s];
    vector<Operation> operations = shift.getOperations();
    
    
    for(int o=0; o<operations.size(); o++){
      bool flag = false;
      vector<unsigned int> allowedTrailers = this->customers[operations[o].getPoint()].getAllowedTrailers();
      //Every trailer is allowed at the source location
      if(operations[o].getPoint()==1)
        flag=true;
      for(int t=0; t<allowedTrailers.size(); t++ ){
	
        if( shift.getTrailer()==allowedTrailers[t] )
          flag = true;
      }
      if(not(flag)){
	 cout<<"SHI05 NOT FEASIBLE!!!!";
	 /* exit(0);*/
	  return true;
      }
    }
   }

  return false;   
}

/**
 * Shift Constraint:
 * [SHI06 | trailerQuantity cannot be negative or exceed capacity of the trailer]
 * [SHI07 | trailerQuantity of a trailer for a shift is the end quantity of the trailer following the previous shift]
 */
bool Instance::shi06(irpSolution solution){

  vector< vector<double> > trailerQuantities(this->trailers.size());
  for(int t=0; t<this->trailers.size(); t++){
    vector<double> trailerQuantity;
    trailerQuantity.push_back(this->trailers[t].getInitialQuantity());
    trailerQuantities[t]=trailerQuantity;
  }

  for(int t=0; t<this->trailers.size(); t++){  
    unsigned int oper = 0;
    
    for(int s=0; s<solution.getShifts().size(); s++){
      Shift shift=solution.getShifts()[s];
      vector<Operation> operations = shift.getOperations();

      if(this->trailers[t].getIndex()==shift.getTrailer()){
        for(int o=0; o<operations.size(); o++){
          if(oper==0)
            trailerQuantities[shift.getTrailer()][oper] -= operations[o].getQuantity();
          else
            trailerQuantities[shift.getTrailer()].push_back(trailerQuantities[shift.getTrailer()][oper-1] - operations[o].getQuantity());

          if( not(
                      trailerQuantities[shift.getTrailer()][oper]>=-EPSILON
                  and trailerQuantities[shift.getTrailer()][oper]<=this->trailers[shift.getTrailer()].getCapacity() + EPSILON   )
                  ){
            cout<<"SHI06 NOT FEASIBLE!!!!  ";
             cout<<s<<" "<<o<<"     "<<oper<<"    "<<trailerQuantities[shift.getTrailer()][oper]<<"\n";
           /* exit(0);*/
            return true;
          }
          oper++;

        }
      }
    }
  } 
  
  return false;  
}

double Instance::checkFeasibility(irpSolution solution, double feasibility){
 
    //replace OR with AND
/*  return(
      this->dri01(solution) ||
      this->dri03(solution) ||
      this->dri08(solution) ||
      this->tl01(solution)  ||
      this->tl03(solution)  ||
      this->dyn01(solution, false) ||
      this->shi02(solution) ||
      this->shi05(solution) ||
      this->shi06(solution)
              );
              */
    double feas = this->dyn01(solution, false);
    if(
          this->dri01(solution) ||
          this->dri03(solution) ||
          this->dri08(solution) ||
          this->tl01(solution)  ||
          this->tl03(solution)  ||
          feas ||
          this->shi02(solution) ||
          this->shi05(solution) ||
          this->shi06(solution)
                  )
        return feas;
    else
        return 0;
}

double Instance::computeObjective(irpSolution solution){
  
  double totalQuantity = 0.0;
  double LR = 0.0;
  for(int s=0; s<solution.getShifts().size(); s++){
    
    double travelDist = 0.0;
    double cost = 0.0;
    Shift shift = solution.getShifts()[s];
    vector<Operation> operations = shift.getOperations();
    
    travelDist += this->distMatrices[0][operations.front().getPoint()];

    for(int o=0; o<operations.size(); o++){
      
      if(operations[o].getQuantity()>0)
        totalQuantity += operations[o].getQuantity();
      
      if(o>0)
        travelDist += this->distMatrices[operations[o-1].getPoint()][operations[o].getPoint()];
    }
    travelDist+=this->distMatrices[operations.back().getPoint()][0];
   
    double end = operations.back().getArrival() + this->customers[operations.back().getPoint()].getSetupTime() + this->timeMatrices[operations.back().getPoint()][0];
    cost = this->trailers[shift.getTrailer()].getDistanceCost()*travelDist + this->drivers[shift.getDriver()].getTimeCost()*(end - shift.getStart());
    LR+=cost;
  }
  LR/=totalQuantity;

  return LR;
}

void Instance::sort(vector<unsigned int> &priorities, vector<double> urgency){
  unsigned int app;
  double temp;
    for (int i=0; i<priorities.size(); i++){
        for (int j=i+1; j<priorities.size(); j++){

            if (urgency[i] > urgency[j]){
                app =  priorities[i];
                priorities[i] = priorities[j];
                priorities[j] = app;
		
                temp =  urgency[i];
                urgency[i] = urgency[j];
                urgency[j] = temp;
            }
        }
    }
}


void Instance::sortPair(vector< pair< pair<unsigned int,unsigned int>, unsigned int > > &priorities){
  pair< pair<unsigned int,unsigned int>, unsigned int > app;
    for (int i=0; i<priorities.size(); i++){
        for (int j=i+1; j<priorities.size(); j++){

            if (priorities[i].first.first > priorities[j].first.first){
                app =  priorities[i];
                priorities[i] = priorities[j];
                priorities[j] = app;
            }
        }
    }
}
/*
void quickSort(vector<unsigned int> &priorities, vector<double> &urgency, int left, int right) {

      int i = left, j = right;
      int tmp, tmp2;
      double pivot = urgency[(left + right) / 2];

      while (i <= j) {
            while (urgency[i] < pivot)
                  i++;
            while (urgency[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = priorities[i];
                  priorities[i] = priorities[j];
                  priorities[j] = tmp;
                  tmp = urgency[i];
                  urgency[i] = urgency[j];
                  urgency[j] = tmp;
                  i++;
                  j--;
            }
      };


      if (left < j)
            quickSort(priorities, urgency, left, j);
      if (i < right)
            quickSort(priorities, urgency, i, right);

}
*/

vector<Operation> Instance::recursiveRandomShift(irpSolution solution,
                                                 unsigned int time,
                                                 unsigned int cumTime,
                                                 vector< vector<double> > &horizonQuantities,
                                                 vector<double> &tankQuantities,
                                                 vector<double> &trailerQuantities,
                                                 vector<double> &pastQuantities,
                                                 vector<unsigned int>  &lastOperations,
                                                 vector< pair<pair<unsigned int,unsigned int>,unsigned int> > &timeWindows,
                                                 vector<unsigned int> &customerList,
                                                 vector<Operation> &operations/*,
                                                 double timeWeight,
                                                 double quantityWeight*/,
                                                 vector<double> &deliveredQuantities,
                                                 vector<double> totalForecasts,
                                                 double timeWeight,
                                                 double quantityWeight,
                                                 int ties,
                                                 vector<bool> refuelFlags){

  if(customerList.size() <= 1){
    return operations;
  }
  else{
   unsigned int prec, succ;
   prec = 0;
   succ = 1;

  unsigned int twD=0;
  bool lastDShift = true;
  for(int t=1; t<timeWindows.size(); t++)
    if(timeWindows[t].second == timeWindows[0].second){
      twD = t;
      lastDShift = false;
      break;
    }
  unsigned int twT=0;
  bool lastTShift = true;
  for(int t=1; t<timeWindows.size(); t++)
    if(this->drivers[timeWindows[t].second].getTrailer() == this->drivers[timeWindows[0].second].getTrailer()){
      twT = t;
      lastTShift = false;
      break;
    }
    bool allowedTrailerFlag = false;
    for(int t=0; t<this->customers[customerList[succ]].getAllowedTrailers().size(); t++)
      if(this->customers[customerList[succ]].getAllowedTrailers()[t] == this->drivers[timeWindows[0].second].getTrailer()  or customerList[succ] <= 1)
    allowedTrailerFlag = true;

    if(trailerQuantities[this->drivers[timeWindows[0].second].getTrailer()] <= EPSILON and
            not(refuelFlags[this->drivers[timeWindows[0].second].getTrailer()]))
              customerList.insert(customerList.begin() + 1, 1);


    /////
    double pastQuantity = 0.0;
    double horizon = 0.0;
        if(customerList[succ]!=1){
            pastQuantity = pastQuantities[customerList[succ]];
            pastQuantity += this->customers[customerList[succ]].getInitialTankQuantity();
            pastQuantity += horizonQuantities[customerList[succ]][(int)(time + this->timeMatrices[customerList[prec]][customerList[succ]])/60];
            tankQuantities[customerList[succ]] = pastQuantity;
            if(tankQuantities[customerList[succ]] < 0.0)
                tankQuantities[customerList[succ]] = 0.0;
            horizon = 0.0;
            if(customerList[succ] != 1){
               horizon = pastQuantity + horizonQuantities[customerList[succ]].back() - horizonQuantities[customerList[succ]][(int)(time + this->timeMatrices[customerList[prec]][customerList[succ]])/60];
               horizon -= this->customers[customerList[succ]].getSafetyLevel();
            }
        }
    //////

   while(not(
     cumTime + this->timeMatrices[customerList[prec]][customerList[succ]] + this->customers[customerList[succ]].getSetupTime() + this->timeMatrices[customerList[succ]][0] <= this->drivers[timeWindows[0].second].getMaxDrivingDuration()
     and timeWindows[0].first.first + cumTime + this->timeMatrices[customerList[prec]][customerList[succ]] + this->customers[customerList[succ]].getSetupTime() + this->timeMatrices[customerList[succ]][0] <= timeWindows[0].first.second
     and (timeWindows[twD].first.first - (cumTime + timeWindows[0].first.first + this->timeMatrices[customerList[prec]][customerList[succ]] + this->customers[customerList[succ]].getSetupTime() + this->timeMatrices[customerList[succ]][0]) >= drivers[timeWindows[0].second].getMinInterSHIFTDURATION()
       or lastDShift   )
     and (timeWindows[twT].first.second <= timeWindows[0].first.first or timeWindows[twT].first.first >= timeWindows[0].first.second
      or lastTShift    )
     and (allowedTrailerFlag or customerList[succ]==1)
     )
     or (tankQuantities[customerList[succ]] >= this->customers[customerList[succ]].getCapacity()- EPSILON and customerList[succ] != 1)
     or (lastOperations[customerList[succ]] >= time + this->timeMatrices[customerList[prec]][customerList[succ]]  and customerList[succ] != 1)
     or (horizon >= -EPSILON and customerList[succ] != 1)
  ){

     if(customerList[succ] == 1)
         return operations;
     customerList.erase(customerList.begin() + succ);

     if(customerList.size() <=1 ){
       return operations;
     }

    allowedTrailerFlag = false;
    for(int t=0; t<this->customers[customerList[succ]].getAllowedTrailers().size(); t++)
      if(this->customers[customerList[succ]].getAllowedTrailers()[t] == this->drivers[timeWindows[0].second].getTrailer() or customerList[succ] <= 1)
    allowedTrailerFlag = true;

    //////
        if(customerList[succ]!=1){
            pastQuantity = pastQuantities[customerList[succ]];
            pastQuantity += this->customers[customerList[succ]].getInitialTankQuantity();
            pastQuantity += horizonQuantities[customerList[succ]][(int)(time + this->timeMatrices[customerList[prec]][customerList[succ]])/60];
            tankQuantities[customerList[succ]] = pastQuantity;
            if(tankQuantities[customerList[succ]] < 0.0)
                tankQuantities[customerList[succ]] = 0.0;
            horizon = 0.0;
            if(customerList[succ] != 1){
               horizon = pastQuantity + horizonQuantities[customerList[succ]].back() - horizonQuantities[customerList[succ]][(int)(time + this->timeMatrices[customerList[prec]][customerList[succ]])/60];
               horizon -= this->customers[customerList[succ]].getSafetyLevel();
             }
         }
    //////
    }
    unsigned int driver = timeWindows[0].second;
    unsigned int trailer = this->drivers[driver].getTrailer();
    unsigned int customer = customerList[succ];

    Operation operation;
    operation.setArrival(time + this->timeMatrices[customerList[prec]][customer]);
    /////
    if(lastOperations[customer] < time + this->timeMatrices[customerList[prec]][customer])
        lastOperations[customer] = time + this->timeMatrices[customerList[prec]][customer];
    /////
    operation.setPoint(customer);

    bool refuelFlag = false;

    cout<<"         Tank: "<<tankQuantities[customer]<<"   trail:"<<trailerQuantities[trailer]<<"   horizon: "<<horizon<<" "<<this->customers[customer].getCapacity()<< "\n";

    if(customer == 1){
      operation.setQuantity(-(this->trailers[trailer].getCapacity() - trailerQuantities[trailer]));
      pastQuantities[customer] += operation.getQuantity();
      trailerQuantities[trailer] = - operation.getQuantity();
      refuelFlags[trailer] = false;
    }
    else if(trailerQuantities[trailer] <= fabs(horizon)){
      double quantity = 0.0;
      if(trailerQuantities[trailer] <= this->customers[customer].getCapacity() - tankQuantities[customer]){
          quantity = trailerQuantities[trailer];
          refuelFlags[trailer] = true;
      }
      else
          quantity = this->customers[customer].getCapacity() - tankQuantities[customer];


      pastQuantities[customer] += quantity;
      operation.setQuantity(quantity);
      trailerQuantities[trailer] -= operation.getQuantity();

    }
    else{
      double quantity = 0.0;
      if(fabs(horizon) < this->customers[customer].getCapacity() - tankQuantities[customer])
          quantity = fabs(horizon);
      else
          quantity = this->customers[customer].getCapacity() - tankQuantities[customer];

      pastQuantities[customer] += quantity;
      operation.setQuantity(quantity);
      trailerQuantities[trailer] -= operation.getQuantity();
    }

    unsigned int currentCumTime = this->timeMatrices[customerList[prec]][customer] + this->customers[customer].getSetupTime();
    unsigned int currentSetupTime = this->customers[customer].getSetupTime();

    customerList.erase(customerList.begin());
    customerList.erase(customerList.begin());

    ////////////////////
    vector<unsigned int> customerList;
    for(int c=0; c<this->customers.size(); c++)
    customerList.push_back(c);

    vector<double> urgency(this->customers.size(),0.0);

    for(int c=2; c<this->customers.size(); c++){
      Customer customer = this->customers[c];
      vector<double> forecast = customer.getForecast();
      double totForecast = 0;
      bool flag = true;

      double deliveredQuantity = deliveredQuantities[c];

      unsigned int currentPoint = 0;
      if(solution.getShifts().size() > 0)
          if(solution.getShifts().back().getOperations().size() > 0)
            currentPoint = solution.getShifts().back().getOperations().back().getPoint();

      double totalDistance = 0.0;
      for(int d=0; d<this->customers.size(); d++)
          totalDistance += this->timeMatrices[currentPoint][d];

      double totalForecast = totalForecasts[c];

      unsigned int f = 0;
      double timeFactor, quantityFactor, distanceFactor;
      while(f < forecast.size() and flag){
        totForecast += forecast[f];
        if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
          timeFactor = ((double)f/forecast.size()) * timeWeight/2;
          quantityFactor = (abs(totalForecast - deliveredQuantity)/totalForecast) * quantityWeight/2;
          distanceFactor = 100.0 * EPSILON * (double)this->timeMatrices[currentPoint][c]/(totalDistance);
          distanceFactor *= ties;
          urgency[c] = 1.0 + timeFactor + quantityFactor + distanceFactor;
          flag = false;
  //        cout<<"Fs: "<<(double)this->timeMatrices[currentPoint][c]<<" "<<totalDistance<<"\n";
  //        cout<<"Fs: "<<c<<" "<<timeFactor<<" "<<quantityFactor<<" "<<distanceFactor<<" "<<urgency[c]<<"\n\n";
//            int a;cin>>a;
        }
        f++;
      }
    }

    for(int u=0; u<urgency.size(); u++){
      if(urgency[u] <= 0){
        customerList.erase(customerList.begin() + u);
        urgency.erase(urgency.begin() + u);
        u--;
      }
      if(customerList[u]==customer){
          customerList.erase(customerList.begin() + u);
          urgency.erase(urgency.begin() + u);
          u--;
      }
    }

    this->sort(customerList, urgency);


    ///////////
    customerList.insert(customerList.begin(), customer);

    operations.push_back(operation);
    deliveredQuantities[customer] += operation.getQuantity();
    if(refuelFlags[trailer]){
      customerList.insert(customerList.begin() + 1, 1);
    }

    cout<<"     Operation: "<<operation.getPoint()<<"   arr: "<<operation.getArrival()<<"   q: "<<operation.getQuantity()<<"\n";
    /*
    for(int c=0; c<customerList.size(); c++)
      cout<<customerList[c]<<" ";
    cout<<"\n";

    for(int u=0; u<urgency.size(); u++)
      cout<<urgency[u]<<" ";
    cout<<"\n";
    */
//            int a;cin>>a;
    return recursiveRandomShift(
      solution,
      operation.getArrival() + currentSetupTime,
      cumTime + currentCumTime,
      horizonQuantities, tankQuantities, trailerQuantities, pastQuantities, lastOperations,
      timeWindows, customerList, operations,
      deliveredQuantities,
      totalForecasts,
      timeWeight,
      quantityWeight,
      ties,
      refuelFlags);
  }

}


irpSolution Instance::recursiveRandomSolution(irpSolution solution,
                                           unsigned int time,
                                           vector< vector<double> > &horizonQuantities,
                                           vector<double> &tankQuantities,
                                           vector<double> &trailerQuantities,
                                           vector<double>  &pastQuantities,
                                           vector<unsigned int>  &lastOperations,
                                           vector< pair<pair<unsigned int,unsigned int>,unsigned int> > &timeWindows,
                                           double timeWeight,
                                           double quantityWeight,
                                           int ties,
                                           vector<double> &deliveredQuantities,
                                           double totalDistance, vector<double> &totalForecasts){

  Shift shift;
  shift.setIndex(solution.getShifts().size());
      
    if(timeWindows.size() == 0){
        return solution;
    }
    else{

      shift.setDriver(timeWindows[0].second);
      shift.setTrailer(this->drivers[timeWindows[0].second].getTrailer());
      shift.setStart(timeWindows[0].first.first);
     
      time += timeWindows[0].first.first;
          
      vector<unsigned int> customerList;      
      for(int c=0; c<this->customers.size(); c++)
	  customerList.push_back(c);
      vector<double> urgency(this->customers.size(),0.0);


      for(int c=2; c<this->customers.size(); c++){
        Customer customer = this->customers[c];
        vector<double> forecast = customer.getForecast();
        double totForecast = 0;
        bool flag = true;

        double deliveredQuantity = deliveredQuantities[c];
        unsigned int currentPoint = 0;

        totalDistance = 0.0;
        for(int d=0; d<this->customers.size(); d++)
            totalDistance += this->timeMatrices[currentPoint][d];

        double totalForecast = totalForecasts[c];

        unsigned int f = 0;
        double timeFactor, quantityFactor, distanceFactor;
        while(f < forecast.size() and flag){
          totForecast += forecast[f];
          if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
            timeFactor = ((double)f/forecast.size()) * timeWeight/2;
            quantityFactor = (abs(totalForecast - deliveredQuantity)/totalForecast) * quantityWeight/2;
            distanceFactor = 100.0 * EPSILON * (double)this->timeMatrices[currentPoint][c]/(totalDistance);
            distanceFactor *= ties;
            urgency[c] = 1.0 + timeFactor + quantityFactor + distanceFactor;
            flag = false;
//            cout<<"Fs: "<<(double)this->timeMatrices[currentPoint][c]<<" "<<totalDistance<<"\n";
//            cout<<"Fs: "<<c<<" "<<timeFactor<<" "<<quantityFactor<<" "<<distanceFactor<<" "<<urgency[c]<<"\n\n";
//            int a;cin>>a;
          }
          f++;
        }
      }
	
      for(int u=0; u<urgency.size(); u++)
        if(urgency[u] <= 0){
          customerList.erase(customerList.begin() + u);
          urgency.erase(urgency.begin() + u);
          u--;
        }
	  
      this->sort(customerList, urgency);
//     quickSort(customerList, urgency, 0, customerList.size());
      customerList.insert(customerList.begin(), 0);

/*
      for(int c=0; c<customerList.size(); c++)
        cout<<customerList[c]<<" ";
      cout<<"\n";
	
      for(int u=0; u<urgency.size(); u++)
        cout<<urgency[u]<<" ";
      cout<<"\n";  
*/
//      int a;cin>>a;
      cout<<"SHIFT: "<<shift.getIndex()<<" start: "<<shift.getStart()<<"  \n";
      vector<bool> refuelFlags(this->trailers.size(), false);
      vector<Operation> operations;
      operations = recursiveRandomShift(solution, shift.getStart(), 0, horizonQuantities, tankQuantities, trailerQuantities,
                                        pastQuantities, lastOperations, timeWindows, customerList, operations,
                                        deliveredQuantities,
                                        totalForecasts,
                                        timeWeight,
                                        quantityWeight,
                                        ties,
                                        refuelFlags);
      if(operations.size() > 0 ){
         shift.setOperations(operations);
         vector <Shift> shifts = solution.getShifts();
         shifts.push_back(shift);
         solution.setShifts(shifts);
      }
      
      timeWindows.erase(timeWindows.begin());
      return recursiveRandomSolution(solution, time, horizonQuantities, tankQuantities, trailerQuantities,
                                     pastQuantities, lastOperations, timeWindows,
                                     timeWeight, quantityWeight, ties,
                                     deliveredQuantities, totalDistance, totalForecasts);
    }
}

irpSolution Instance::backTrackingRandomSolution(double timeWeight, double quantityWeight, int ties){
  
  irpSolution solution;
  vector< vector<double> > horizonQuantities(this->customers.size());
  for(int p=0; p<this->customers.size();p++){
    if(p>1){
      vector<double> horizon(this->customers[p].getForecast().size(),0.0);
      horizonQuantities[p]=horizon;
      horizonQuantities[p][0] = 0.0 /*this->customers[p].getInitialTankQuantity()*/;
      for(int f=0; f<this->customers[p].getForecast().size(); f++){
        if(f>0){
          horizonQuantities[p][f] = horizonQuantities[p][f] + horizonQuantities[p][f-1] - this->customers[p].getForecast()[f];
        }
        else{
          horizonQuantities[p][f] -= this->customers[p].getForecast()[f];
        }
      }

    }
    else{
      vector<double> horizon(this->customers[p].getForecast().size(),0.0);
      horizonQuantities[p]=horizon;
    }
  }

  vector<double> tankQuantities(this->customers.size());
  for(int c=0; c<this->customers.size(); c++){
    if(c == 0)
       tankQuantities[c] = 0.0;
    else
      tankQuantities[c] = this->customers[c].getInitialTankQuantity();
  }
  
  vector<double> trailerQuantities(this->trailers.size());
  for(int t=0; t<this->trailers.size(); t++){
    trailerQuantities[t] = this->trailers[t].getInitialQuantity();
  }
  
   vector< pair< pair<unsigned int,unsigned int>, unsigned int > > timeWindows;
   for(int d=0; d<this->drivers.size(); d++){
     for(int t=0; t<this->drivers[d].getTimeWindows().size(); t++){
       pair< pair<unsigned int,unsigned int>, unsigned int > tw;
       tw.first = this->drivers[d].getTimeWindows()[t];
       tw.second = d;
       timeWindows.push_back(tw);
     }
   }
  this->sortPair(timeWindows);
 
  vector<unsigned int> lastOperations(this->customers.size(), 0.0);
  vector<double> pastQuantities(this->customers.size(), 0.0);

  vector<double> deliveredQuantities(this->customers.size(), 0.0);

  /////
  double totalDistance = 0.0;
  for(int d=0; d<this->customers.size(); d++)
      totalDistance += this->timeMatrices[0][d];
  /////

  vector<double> totalForecasts(this->customers.size(), 0.0);
  for(int c=0; c<totalForecasts.size(); c++)
      for(int f=0; f<this->customers[c].getForecast().size(); f++)
          totalForecasts[c] += this->customers[c].getForecast()[f];

  
  return recursiveRandomSolution(solution,
                                 0,
                                 horizonQuantities, tankQuantities, trailerQuantities,
                                 pastQuantities, lastOperations,
                                 timeWindows, timeWeight, quantityWeight, ties,
                                 deliveredQuantities,
                                 totalDistance,
                                 totalForecasts);
}



irpSolution Instance::rebuildSolution(irpSolution &initialSolution, vector<unsigned int> solutionRepresentation,
                                      double refuelRatio, double deliveredQuantityRatio){

    deliveredQuantityRatio = 1.0 - deliveredQuantityRatio;
/*    cout<<"REPRESENTATION: ";
    for(int p=0; p<solutionRepresentation.size(); p++)
        cout<<solutionRepresentation[p]<<" ";
    cout<<"\n\n";
*/
    vector< vector<double> > horizons;
    vector< vector<double> > actualQuantity;
    for(int c=0; c<this->customers.size(); c++){
        vector<double> horizon(this->horizon*60, 0.0);
        horizons.push_back(horizon);
        actualQuantity.push_back(horizon);
        horizons[c][0] += this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel();
        actualQuantity[c][0] += this->customers[c].getInitialTankQuantity();
        if(c>1)
            for(int f=0; f<this->horizon*60; f++){
                if(f>0){
     /*               horizons[c][f] += (horizons[c][f-1] - this->customers[c].getForecast()[(int)f/60]);
                    actualQuantity[c][f] += (actualQuantity[c][f-1] - this->customers[c].getForecast()[(int)f/60]);*/
                    horizons[c][f] += horizons[c][f-1];
                    actualQuantity[c][f] += actualQuantity[c][f-1];
                    if(f%60 == 0){
                        horizons[c][f] -= this->customers[c].getForecast()[(int)f/60];
                        actualQuantity[c][f] -= this->customers[c].getForecast()[(int)f/60];
                    }
                }
                else{
                    horizons[c][f] -= this->customers[c].getForecast()[0];
                    actualQuantity[c][f] -= this->customers[c].getForecast()[0];
                }
            }
    }

    vector< pair< pair<unsigned int,unsigned int>, unsigned int > > timeWindows;
    for(int d=0; d<this->drivers.size(); d++){
      for(int t=0; t<this->drivers[d].getTimeWindows().size(); t++){
        pair< pair<unsigned int,unsigned int>, unsigned int > tw;
        tw.first = this->drivers[d].getTimeWindows()[t];
        tw.second = d;
        timeWindows.push_back(tw);
      }
    }
   this->sortPair(timeWindows);

   irpSolution solution;
   vector <Shift> shifts;
   unsigned int time;
   unsigned int index = 0;

   vector<double> trailerQuantities(this->trailers.size(), 0.0);
   for(int t=0; t<trailerQuantities.size(); t++)
       trailerQuantities[t] = this->trailers[t].getInitialQuantity();
   vector<double> tankQuantities(this->customers.size(), 0.0);
   for(int c=0; c<tankQuantities.size(); c++)
       tankQuantities[c] = this->customers[c].getInitialTankQuantity();

   vector <bool> refuelFlag(this->trailers.size(), false);
   for(int t=0; t<this->trailers.size(); t++)
       if(trailerQuantities[t] <= EPSILON)
           refuelFlag[t] = true;
   while(solutionRepresentation.size() > 0 and timeWindows.size() > 0){

       unsigned int driver = timeWindows[0].second;
       unsigned int trailer = this->drivers[driver].getTrailer();
       Shift shift;
       shift.setIndex(index);
       shift.setDriver(driver);
       shift.setTrailer(this->drivers[timeWindows[0].second].getTrailer());
       shift.setStart(timeWindows[0].first.first);


       for(int t=0; t<this->trailers.size(); t++)
           if(trailerQuantities[t] <= EPSILON)
               refuelFlag[t] = true;
/*
       cout<<"SHIFT: "<<shift.getIndex()<<"\n";
       cout<<"       "<<shift.getDriver()<<" "<<shift.getTrailer()<<" "<<shift.getStart()<<"\n";
*/
//       cout<<"SHIFT: "<<shift.getIndex()<<" start: "<<shift.getStart()<<" "<<shift.getDriver()<<" "<<shift.getTrailer()<<"\n";
       time = timeWindows[0].first.first;

       unsigned int prec;
       unsigned int succ;
       if(refuelFlag[trailer]){
           prec = 0;
           succ = 1;
       }
       else{
           prec = 0;
           succ = solutionRepresentation.front();
       }

       unsigned int cumulatedTime = time + this->timeMatrices[0][succ];

       unsigned int twD=0;
       bool lastDShift = true;
       for(int t=1; t<timeWindows.size(); t++)
         if(timeWindows[t].second == driver){
           twD = t;
           lastDShift = false;
           break;
         }
       unsigned int twT=0;
       bool lastTShift = true;
       for(int t=1; t<timeWindows.size(); t++)
         if(this->drivers[timeWindows[t].second].getTrailer() == this->drivers[driver].getTrailer()){
           twT = t;
           lastTShift = false;
           break;
         }
         bool allowedTrailerFlag = false;
         for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
           if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[driver].getTrailer())
               allowedTrailerFlag = true;

       vector<Operation> operations;
       bool horizonFlag = false;

       while((cumulatedTime + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0] - time) <= this->drivers[driver].getMaxDrivingDuration()
              and cumulatedTime + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0] <= timeWindows[0].first.second
             and solutionRepresentation.size() > 0
             and (timeWindows[twD].first.first - (cumulatedTime + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) >= this->drivers[driver].getMinInterSHIFTDURATION()
                  or lastDShift)
             and (timeWindows[twT].first.second <= timeWindows[0].first.first or timeWindows[twT].first.first >= timeWindows[0].first.second
              or lastTShift    )
             and (allowedTrailerFlag or succ == 1)
             and timeWindows.size()!=0

           ){
           Operation operation;
           operation.setArrival(cumulatedTime);
           operation.setPoint(succ);

//           cout<<"QU; "<<trailerQuantities[trailer]<<" "<<actualQuantity[succ][operation.getArrival()]/*tankQuantities[succ]*/<<"     "<<horizons[succ].back()<<" "<<this->customers[succ].getCapacity()<<"       ";
 //          cout<<"         Tank: "<<actualQuantity[succ][operation.getArrival()]<<"   trail:"<<trailerQuantities[trailer]/*<<"   horizon: "<<horizons[succ].back()*/<<" "<<this->customers[succ].getCapacity()<<"       ";
           double quantity = 0.0;

           if(refuelFlag[trailer] == true){
               quantity = this->trailers[trailer].getCapacity() - trailerQuantities[trailer];
               trailerQuantities[trailer] += quantity;
               refuelFlag[trailer] = false;

               operation.setQuantity(-quantity);

       //        cout<<"OP: "<<operation.getPoint()<<" "<<operation.getArrival()<<" "<<operation.getQuantity()<<"\n";
  //             cout<<"     Operation: "<<operation.getPoint()<<"   arr: "<<operation.getArrival()<<"   q: "<<operation.getQuantity()<<"\n";
               operations.push_back(operation);

               prec = 1;
               succ = solutionRepresentation.front();

               cumulatedTime += this->customers[prec].getSetupTime() + this->timeMatrices[prec][succ];

               allowedTrailerFlag = false;
               for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                 if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[driver].getTrailer())
                     allowedTrailerFlag = true;
           }
           else{
           if(trailerQuantities[trailer] <= fabs(horizons[succ].back())
                   and horizons[succ].back() < 0){
               double residualQuantity = actualQuantity[succ][(int)(operation.getArrival()/* + this->customers[succ].getSetupTime()*/)/*/60*/];
               if(residualQuantity < 0)
                   residualQuantity = 0;
               else if(residualQuantity > this->customers[succ].getCapacity())
                   residualQuantity = this->customers[succ].getCapacity();


               if(trailerQuantities[trailer] <= this->customers[succ].getCapacity() - residualQuantity/*tankQuantities[succ*/){
                   quantity = trailerQuantities[trailer] /** deliveredQuantityRatio*/;

                   for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                       if(actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                          quantity -= (actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                          horizonFlag = true;
                       }
                   }
                   if(horizonFlag)
                       break;

                   /////
                   double factor = 1.0/* - residualQuantity/this->customers[succ].getCapacity()*/;
                   if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                   quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                   /////

                   trailerQuantities[trailer] -= quantity;
                   tankQuantities[succ] += quantity;
                   /*
                   if(trailerQuantities[trailer] < EPSILON)
                     refuelFlag[trailer] = true;
                     */
               }
               else{
                   quantity = (this->customers[succ].getCapacity() - residualQuantity) /** deliveredQuantityRatio*/;

                   for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                       if(actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                          quantity -= (actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                          horizonFlag = true;
                       }
                   }
                   if(horizonFlag)
                       break;

                   /////
                   double factor = 1.0/* - residualQuantity/this->customers[succ].getCapacity()*/;
                   if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                       quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                   /////

                   trailerQuantities[trailer] -= quantity;
                   tankQuantities[succ] += quantity;
               }
           }
           else{
               double residualQuantity = actualQuantity[succ][(int)(operation.getArrival()/* + this->customers[succ].getSetupTime()*/)/*/60*/];
               if(residualQuantity < 0)
                   residualQuantity = 0;
               else if(residualQuantity > this->customers[succ].getCapacity())
                   residualQuantity = this->customers[succ].getCapacity();

               if(fabs(horizons[succ].back()) <= this->customers[succ].getCapacity() - residualQuantity/*tankQuantities[succ]*/
                       and horizons[succ].back() < 0){
                   quantity = fabs(horizons[succ].back()) /** deliveredQuantityRatio*/;

                   for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                       if(actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                          quantity -= (actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                          horizonFlag = true;
                       }
                   }
                   if(horizonFlag)
                       break;

                   /////
                   double factor = 1.0/* - residualQuantity/this->customers[succ].getCapacity()*/;
                   if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                      quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                   /////

                   trailerQuantities[trailer] -= quantity;
                   tankQuantities[succ] += quantity;
               }
               else if(fabs(horizons[succ].back()) > this->customers[succ].getCapacity() - residualQuantity/*tankQuantities[succ]*/
                       and horizons[succ].back() < 0){
                   quantity = (this->customers[succ].getCapacity() - residualQuantity) /** deliveredQuantityRatio*/;

                   for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                       if(actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                          quantity -= (actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                          horizonFlag = true;
                        }
                   }
                   if(horizonFlag)
                       break;

                   /////
                   double factor = 1.0/* - residualQuantity/this->customers[succ].getCapacity()*/;
                   if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                        quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                   /////

                   trailerQuantities[trailer] -= quantity;
                   tankQuantities[succ] += quantity;

               }
               else
                   quantity = 0.0;
           }

           if(horizonFlag or (quantity <= EPSILON and operation.getPoint()!=1))
               break;
           if(trailerQuantities[trailer] <= EPSILON + this->trailers[trailer].getCapacity() * refuelRatio){
//             cout<<"   REFUEL!!!   ";
             refuelFlag[trailer] = true;
           }

           for(int ff = (int)(operation.getArrival() /*+ this->customers[succ].getSetupTime()*/)/*/60*/; ff < this->customers[succ].getForecast().size()*60; ff++){
               horizons[succ][ff] += quantity;
               actualQuantity[succ][ff] += quantity;
           }

           operation.setQuantity(quantity);


       //    cout<<"OP: "<<operation.getPoint()<<" "<<operation.getArrival()<<" "<<operation.getQuantity()<<"\n";
 //          cout<<"     Operation: "<<operation.getPoint()<<"   arr: "<<operation.getArrival()<<"   q: "<<operation.getQuantity()<<"\n";
           operations.push_back(operation);
           solutionRepresentation.erase(solutionRepresentation.begin());
           if(refuelFlag[trailer]){
 //              solutionRepresentation.insert(solutionRepresentation.begin(), 1);
 //              refuelFlag[trailer] = false;
               prec = succ;
               succ = 1;
           }
           else{
               prec = succ;
               if(solutionRepresentation.size() > 0){
                  /* if(horizons[solutionRepresentation.front()].back() >= 0)
                       solutionRepresentation.erase(solutionRepresentation.begin());*/
                   if(solutionRepresentation.size() > 0)
                      succ = solutionRepresentation.front();
               }
           }
           cumulatedTime += this->customers[prec].getSetupTime() + this->timeMatrices[prec][succ];

           allowedTrailerFlag = false;
           for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
             if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[driver].getTrailer())
                 allowedTrailerFlag = true;
       }

       }
//       cout<<"\n";
       timeWindows.erase(timeWindows.begin());
       shift.setOperations(operations);
       if(operations.size() != 0)
        shifts.push_back(shift);
       else
           index--;

       index++;
   }

   solution.setShifts(shifts);

/*
   cout<<"REBUILT SOLUTION: ";
   for(int s=0; s<solution.getShifts().size();s++){
       cout<<"SHIFT: "<<solution.getShifts()[s].getIndex()<<" "<<solution.getShifts()[s].getStart()<<"\n";
       for(int o=0; o<solution.getShifts()[s].getOperations().size();o++)
           cout<<solution.getShifts()[s].getOperations()[o].getPoint()<<" "<<solution.getShifts()[s].getOperations()[o].getArrival()<<" "<<solution.getShifts()[s].getOperations()[o].getQuantity()<<"\n";
       cout<<"\n";
   }

   cout<<"HORIZON: \n";
   for(int h=0; h<horizons.size(); h++)
       cout<<horizons[h][0]<<"   "<<horizons[h].back()<<"\n";
    int a;
    */
/*
   cout<<"HORIZON: \n";
   for(int f=0; f<this->horizon*60; f++){
       cout<<actualQuantity[6][f]<<"\n";
       if(actualQuantity[6][f] > 6000){
                      cout<<f<<" "<<actualQuantity[6].back()<<"\n";
           cin>>a;
       }
   }
*/

   return solution;
}


void Instance::compute(irpSolution &solution){

    cout<<"INITIAL TRAILER QUANTITY: ";
    for(int t=0; t<this->trailers.size(); t++){
        vector<double> trailerQuantity(this->horizon * 60,0.0);
        solution.trailerQuantities.push_back(trailerQuantity);
        solution.trailerQuantities[t][0] = this->trailers[t].getInitialQuantity();
        cout<<solution.trailerQuantities[t][0]<<" ";
    }
    cout<<"\n";
    cout<<"INITIAL TANK QUANTITY: ";
    for(int c=0; c<this->customers.size(); c++){
        vector<double> tankQuantity(this->horizon * 60,0.0);
        solution.tankQuantities.push_back(tankQuantity);
        if(c > 1)
            solution.tankQuantities[c][0] = this->customers[c].getInitialTankQuantity();
//        cout<<solution.tankQuantities[c][0]<<" ";
    }

    for(int c=0; c<this->customers.size(); c++){
        for(int f=0; f<this->customers[c].getForecast().size(); f++){
            if(c > 1){
                solution.tankQuantities[c][f*60] -= this->customers[c].getForecast()[f];
                if(c==2)
                for(int ff=0; ff<60; ff++);
//                    cout<<solution.tankQuantities[c][f*60 + ff]<<" ";
            }
        }
    }

    for(int s=0; s<solution.getShifts().size(); s++){
        for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++){
            Operation op = solution.getShifts()[s].getOperations()[o];
            if(op.getPoint() > 0){
                unsigned int opTime = op.getArrival() + this->customers[op.getPoint()].getSetupTime();
                unsigned int driver = solution.getShifts()[s].getDriver();
                solution.tankQuantities[op.getPoint()][opTime] += op.getQuantity();
                solution.trailerQuantities[driver][opTime] -= op.getQuantity();
            }
        }
    }

    for(int c=0; c<this->customers.size(); c++){
        for(int f=0; f<this->customers[c].getForecast().size(); f++){
            if(c > 0)
               for(int ff=0; ff<60; ff++)
                   if(f*60 + ff > 0)
                       solution.tankQuantities[c][f*60 + ff] += solution.tankQuantities[c][f*60 + ff - 1];
        /*    if(c==2 and f <15)
            for(int ff=0; ff<60; ff++)
                cout<<f<<" "<<f*60+ff<<" "<<solution.tankQuantities[c][f*60 + ff]<<"\n";*/
            }
    }

    for(int t=0; t<this->trailers.size(); t++){
        for(int f=0; f<this->horizon; f++){
               for(int ff=0; ff<60; ff++)
                   if(f*60 + ff > 0)
                       solution.trailerQuantities[t][f*60 + ff] += solution.trailerQuantities[t][f*60 + ff - 1];

          /*  if(t==0 and f <15)
             for(int ff=0; ff<60; ff++)
                   cout<<f<<" "<<f*60+ff<<" "<<solution.trailerQuantities[t][f*60 + ff]<<"\n";*/
            }
    }



}


irpSolution Instance::twoExchangeNeighborhood(irpSolution solution, int s1, int o1, int s2, int o2){

    irpSolution  neighboringSolution = solution;
    unsigned int point1 = neighboringSolution.getShifts()[s1].getOperations()[o1].getPoint();
    unsigned int point2 = neighboringSolution.getShifts()[s2].getOperations()[o2].getPoint();
    if(  not(  (s1==s2 and o1==o2) or (point1==point2)  )  ){
//      cout<<"\n\nCIAOOO "<<neighboringSolution.getShifts()[s2].getOperations()[o2].getPoint();
        neighboringSolution.setPoint(point2, s1, o1);
        neighboringSolution.setPoint(point1, s2, o2);
//      cout<<"\n\nCIAOOO "<<neighboringSolution.getShifts()[s2].getOperations()[o2].getPoint();
        if(!this->checkFeasibility(neighboringSolution, false)){
            for(int s=0; s<solution.getShifts().size(); s++)
                for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++)
                    cout<<neighboringSolution.getShifts()[s].getOperations()[o].getPoint()<<" ";
        }

    }

    return neighboringSolution;
}



#endif

