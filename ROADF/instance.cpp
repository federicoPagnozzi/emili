#ifndef __INSTANCE_CPP
#define __INSTANCE_CPP

#include "instance.h"
#include "tinyxml.h"
#include "tinystr.h"

//#include "/home/antoniofisk/Desktop/Uni/MasterThesis/Project/emilibase.h"
#include "../emilibase.h"

#define COUT if (0) cout
#define CIN if (0) cin

#define EPSILON 0.000001
#define FEASIBILITY_PENALTY 1

Instance::Instance()
{
}

string Instance::getName()  {return this->name;}
int Instance::getUnit()	{return this->unit;}
int Instance::getHorizon()	{return this->horizon;}
vector<vector<unsigned int> > Instance::getTimeMatrices()	{return this->timeMatrices;}
vector<Driver> Instance::getDrivers()	{return this->drivers;}
vector<Trailer> Instance::getTrailers()	{return this->trailers;}
vector<Customer> Instance::getCustomers()	{return this->customers;}
vector<vector<double> > Instance::getDistMatrices()	{return this->distMatrices;}

vector< vector<double> > Instance::getHorizons()    {return this->horizons;}
vector< vector<double> > Instance::getActualQuantity()  {return this->actualQuantity;}
void Instance::setHorizons(vector< vector<double> > h)  {this->horizons = h;}
void Instance::setActualQuantity(vector< vector<double> > aq)   {this->actualQuantity = aq;}

vector< pair< pair<unsigned int,unsigned int>, unsigned int > > Instance::getTimeWindows()  {return this->timeWindows;}
void Instance::setTimeWindows(vector< pair< pair<unsigned int,unsigned int>, unsigned int > > tw)   {this->timeWindows = tw;}

map< unsigned int , vector< vector<unsigned int> >   > Instance::getContributeMatrixes() {return this->contributeMatrixes;}
void Instance::setContributeMatrixes(map<unsigned int, vector<vector<unsigned int> > > cm) {this->contributeMatrixes = cm;}

map< unsigned int , vector< vector<double> >   > Instance::getContributeMatrixesValues() {return this->contributeMatrixesValues;}
void Instance::setContributeMatrixesValues(map<unsigned int, vector<vector<double> > > cm) {this->contributeMatrixesValues = cm;}

vector< vector< vector <unsigned int> > > Instance::getFastestSubroute()    {return this->fastestSubroute;}
void Instance::setFastestSubroute(vector< vector< vector <unsigned int> > > fs) {this->fastestSubroute = fs;}



/**
 * Read xml instance file
 */
void Instance::loadInstance(const char* pFilename)
{
      this->name = pFilename;
    this->name = this->name.substr(this->name.size()-18, this->name.size()-4);
      int index;
      double d;
      TiXmlDocument doc(pFilename);

      if (!doc.LoadFile()) return;
      
      TiXmlHandle hDoc(&doc);
      TiXmlElement* pInnerElem;
      TiXmlHandle hRoot(0);


      pInnerElem=hDoc.FirstChildElement().Element();
      if (!pInnerElem) return;
//      name=pInnerElem->Value();

      hRoot=TiXmlHandle(pInnerElem);
      
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild().Element();
      
      //Read unit
      istringstream ( pInnerElem->GetText()) >> this->unit;
//    COUT<<pInnerElem->GetText()<<" "<<this->unit<<"\n";
      
      //Read horizon
      pInnerElem=pInnerElem->NextSiblingElement();
      istringstream ( pInnerElem->GetText()) >> this->horizon;
//    COUT<<pInnerElem->Value()<<" "<<this->horizon<<"\n";
      
      //Read time matrix
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Instance" ).FirstChild("timeMatrices").FirstChild().Element();
      int i = 0;
      for( pInnerElem; pInnerElem; pInnerElem=pInnerElem->NextSiblingElement())
      {
	      int j = 0;
          vector<unsigned int> vec;
	      this->timeMatrices.push_back(vec);
	      TiXmlElement *pMatrixElem = pInnerElem->FirstChild()->ToElement();
	      for( pMatrixElem ; pMatrixElem;pMatrixElem =pMatrixElem ->NextSiblingElement())
	      {
		int elem;
		istringstream(pMatrixElem ->GetText()) >> elem;
		this->timeMatrices.at(i).push_back(elem);
//		COUT<<this->timeMatrices.at(i).at(j)<<" ";
		j++;
	      }
//	      COUT<<"\n";
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
//	COUT<<pDriverElem->Value()<< " " << d.getIndex()<<"\n";
	
	//Read driver's masimum driving duration
	pDriverElem=pDriverElem->NextSiblingElement();	
	istringstream(pDriverElem->GetText()) >> index;
	d.setMaxDrivingDuration(index);	
//	COUT<<pDriverElem->Value()<< " " << d.getMaxDrivingDuration()<<"\n";
	
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
//	COUT<<pDriverElem->Value()<< " " << d.getTrailer()<<"\n";
	
	//Read driver's minimum interval shift duration
	pDriverElem=pDriverElem->NextSiblingElement();	
	istringstream(pDriverElem->GetText()) >> index;
	d.setMinInterSHIFTDURATION(index);	
//	COUT<<pDriverElem->Value()<< " " << d.getMinInterSHIFTDURATION()<<"\n";
	
	//Read driver's time cost
	double tc;
	pDriverElem=pDriverElem->NextSiblingElement();	
	istringstream(pDriverElem->GetText()) >> tc;
	d.setTimeCost(tc);	
//	COUT<<pDriverElem->Value()<< " " <<d.getTimeCost()<<"\n";
	
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
//	COUT<<pTrailerElem->Value()<< " " << t.getIndex()<<"\n";
	
	//Read trailer's capacity
	pTrailerElem=pTrailerElem->NextSiblingElement();	
	istringstream(pTrailerElem->GetText()) >> d;
	t.setCapacity(d);	
//	COUT<<pTrailerElem->Value()<< " " << t.getCapacity()<<"\n";
	
	//Read trailer's initial quantity
	pTrailerElem=pTrailerElem->NextSiblingElement();	
	istringstream(pTrailerElem->GetText()) >> d;
	t.setInitialQuantity(d);	
//	COUT<<pTrailerElem->Value()<< " " << t.getInitialQuantity()<<"\n";
	
	//Read trailer's distance cost
	pTrailerElem=pTrailerElem->NextSiblingElement();	
	istringstream(pTrailerElem->GetText()) >> d;
	t.setDistanceCost(d);	
//	COUT<<pTrailerElem->Value()<< " " << t.getDistanceCost()<<"\n";
	
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
//	COUT<<pInnerElem->Value()<< " " << b.getIndex()<<"\n";

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
//	COUT<<pSourceElem->Value()<< " " << s.getIndex()<<"\n";
	
	//Read source's setup time
	pSourceElem=pSourceElem->NextSiblingElement();	
	istringstream(pSourceElem->GetText()) >> index;
	s.setSetupTime(index);	
//	COUT<<pSourceElem->Value()<< " " << s.getSetupTime()<<"\n";
	
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
//	COUT<<pCustomerElem->Value()<< " " << c.getIndex()<<"\n";
	
	//Read customer's setup time
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> index;
	c.setSetupTime(index);	
//	COUT<<pCustomerElem->Value()<< " " << c.getSetupTime()<<"\n";
	
        //Read customer's allowed trailers
	pCustomerElem=pCustomerElem->NextSiblingElement();
	TiXmlElement *pAllowedTrailersElem=pCustomerElem->FirstChild()->ToElement();
	vector<unsigned int> at;
	while(pAllowedTrailersElem!= NULL and strcmp(pAllowedTrailersElem->Value(), "int")==0){
	    
	    istringstream(pAllowedTrailersElem->ToElement()->GetText()) >> index;
	    at.push_back(index);
//	    COUT<<pAllowedTrailersElem->Value()<< " " << at.at(0)<<"\n";
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
//	    COUT<<"here here"<<pForecastElem->Value()<< " " << f.at(0)<<"\n";
	    pForecastElem=pForecastElem->NextSiblingElement();
	}
	c.setForecast(f);
	
	//Read customer's capacity
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> d;
	c.setCapacity(d);	
//	COUT<<pCustomerElem->Value()<< " " << c.getCapacity()<<"\n";
	
	//Read customer's initial tank quantity
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> d;
	c.setInitialTankQuantity(d);	
//	COUT<<pCustomerElem->Value()<< " " << c.getInitialTankQuantity()<<"\n";
	
	//Read customer's safety level
	pCustomerElem=pCustomerElem->NextSiblingElement();	
	istringstream(pCustomerElem->GetText()) >> d;
	c.setSafetyLevel(d);	
//	COUT<<pCustomerElem->Value()<< " " << c.getSafetyLevel()<<"\n";
	 
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
//	      COUT<<pInnerElem ->Value()<<": ";
	      TiXmlElement *pMatrixElem = pInnerElem->FirstChild()->ToElement();
	      for( pMatrixElem ; pMatrixElem;pMatrixElem =pMatrixElem ->NextSiblingElement())
	      {
		int elem;
		istringstream(pMatrixElem ->GetText()) >> elem;
		this->distMatrices.at(i).push_back(elem);
    //	COUT<<this->distMatrices.at(i).at(j)<<" ";
		j++;
	      }
    //      COUT<<"\n";
	      i++;
      }


      for(int c=0; c<this->customers.size(); c++){
          vector<double> horizon(this->horizon*60, 0.0);
          this->horizons.push_back(horizon);
          this->actualQuantity.push_back(horizon);
          this->horizons[c][0] += this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel();
          this->actualQuantity[c][0] += this->customers[c].getInitialTankQuantity();
          if(c>1)
              for(int f=0; f<this->horizon*60; f++){
                  if(f>0){
                      this->horizons[c][f] += this->horizons[c][f-1];
                      this->actualQuantity[c][f] += this->actualQuantity[c][f-1];
                      if(f%60 == 0){
                          this->horizons[c][f] -= this->customers[c].getForecast()[(int)f/60];
                          this->actualQuantity[c][f] -= this->customers[c].getForecast()[(int)f/60];
                      }
                  }
                  else{
                      this->horizons[c][f] -= this->customers[c].getForecast()[0];
                      this->actualQuantity[c][f] -= this->customers[c].getForecast()[0];
                  }
              }
      }

//      vector< pair< pair<unsigned int,unsigned int>, unsigned int > > timeWindows;
      for(int d=0; d<this->drivers.size(); d++){
        for(int t=0; t<this->drivers[d].getTimeWindows().size(); t++){
          pair< pair<unsigned int,unsigned int>, unsigned int > tw;
          tw.first = this->drivers[d].getTimeWindows()[t];
          tw.second = d;
          this->timeWindows.push_back(tw);
        }
      }
     this->sortPair(this->timeWindows);


      for(int d=0; d<this->drivers.size(); d++){

          unsigned int driver = this->drivers[d].getIndex();
          unsigned int trailer = this->drivers[driver].getTrailer();
          vector< vector<unsigned int> > closestCustomers(this->customers.size());
          vector< vector<double> > closestCustomersValues(this->customers.size());
          vector<unsigned int> customerList;
          for(int c=0; c<this->customers.size(); c++)
              customerList.push_back(c);
          for(int c1=0; c1<closestCustomers.size(); c1++){
              vector<double> cost(this->customers.size());
              for(int c2=0; c2<cost.size(); c2++){
                  double distanceCost = this->trailers[trailer].getDistanceCost() * this->distMatrices[c1][c2];
                  double timeCost = this->drivers[driver].getTimeCost() * this->timeMatrices[c1][c2];

                  if(c1 == c2)
                      cost[c2] = DBL_MAX;
                  else
                    cost[c2] =  (distanceCost + timeCost);


              }
              closestCustomersValues[c1] = cost;
              this->sort(customerList, cost);
              closestCustomers[c1] = customerList;


          }
          pair< unsigned int, vector< vector<unsigned int> >> p;
          p.first = driver;
          p.second = closestCustomers;
          this->contributeMatrixes.insert(p);

          pair< unsigned int, vector< vector<double> >> pv;
          pv.first = driver;
          pv.second = closestCustomersValues;
          this->contributeMatrixesValues.insert(pv);
      }

      for(int d=0; d<this->drivers.size(); d++){
          COUT<<"DRIVER: "<<this->drivers[d].getIndex()<<"\n";
          vector <vector <unsigned int> > cc = this->contributeMatrixes.find(d)->second;
          for(int c1=0; c1<cc.size(); c1++){
              vector<unsigned int> ccc = cc[c1];
              for(int c2=0; c2<ccc.size(); c2++)
                  COUT<<ccc[c2]<<" ";
              COUT<<"\n";
          }
      }

      vector<vector <unsigned int>> route(this->customers.size(), vector<unsigned int>());

      vector< vector< vector <unsigned int> > > fastestSubroute;
      for(int c1 = 0; c1<this->customers.size(); c1++){
          vector< vector<unsigned int> > v1;
          for(int c2 = 0; c2<this->customers.size(); c2++){
              vector<unsigned int> v2;
              for(int c3 = 0; c3<this->customers.size(); c3++){
                unsigned int time = this->timeMatrices[c1][c3] + this->customers[c3].getSetupTime() + this->timeMatrices[c3][c2];
                v2.push_back(time);
              }
              v1.push_back(v2);
          }
          fastestSubroute.push_back(v1);
      }
      this->fastestSubroute = fastestSubroute;

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
//	  COUT<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<shift1.getStart()<<" "<<end1<<" \n";
	  if(  not(  shift2.getStart() >= end1  + driver.getMinInterSHIFTDURATION() 
          //  or shift1.getStart()  >= end2 + driver.getMinInterSHIFTDURATION()
                 )
           ){
        COUT<<"dri01 NOT FEASIBLE!!!!\n";
//	    COUT<<shift1.getStart()<<" "<<end1<<" "<<shift2.getStart()<<" "<<end2<<"\n";
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
    //	  COUT<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<shift1.getStart()<<" "<<end1<<" \n";
          if(  not(  shift2.getStart() >= end1  + driver.getMinInterSHIFTDURATION()
                /*or shift1.getStart()  >= end2 + driver.getMinInterSHIFTDURATION() */)
               ){
            COUT<<"dri01 NOT FEASIBLE!!!!\n";
            COUT<<s1<<" "<<s2<<"\n";
            COUT<<shift1.getStart()<<" "<<end1<<" "<<shift2.getStart()<<" "<<end2<<"\n";
            int a; CIN>>a;
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
    
// COUT<<"CHECK SHIFT: "<<  solution.getShifts().at(s).getIndex()<<" ";
   vector<Operation> operations = solution.getShifts()[s].getOperations();
   cumulatedDrivingTime += this->timeMatrices[0][operations.front().getPoint()];
   for(int o=1; o<operations.size(); o++){
     
     cumulatedDrivingTime+=this->timeMatrices[operations[o-1].getPoint()][operations[o].getPoint()];
 //    COUT<<operations.at(o-1).getPoint()<<" "<<operations.at(o).getPoint()<<" "<<cumulatedDrivingTime<<"\n";
   }

   cumulatedDrivingTime += this->timeMatrices[operations.back().getPoint()][0];
   if(not(cumulatedDrivingTime <= this->drivers.at(solution.getShifts()[s].getDriver()).getMaxDrivingDuration())){
     COUT<<"\nDRI03 NOT FEASIBLE!!!!";
     COUT<<s<<" "<<cumulatedDrivingTime<<"\n";
     int a; CIN>>a;
     /* exit(0);*/
      return true;
   }
//    COUT<<"\n";
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
//    COUT<<"SHIFT: "<<shift.getIndex()<<"\n";
    for(int tw=0; tw<timeWindows.size(); tw++){
      
      unsigned int end = shift.getOperations().back().getArrival() + this->customers.at(shift.getOperations().back().getPoint()).getSetupTime()
			 + this->timeMatrices[shift.getOperations().back().getPoint()][0];
//      COUT<<timeWindows[tw].first<<" "<<timeWindows[tw].second<<" "<<shift.getStart()<<" "<<end<<"\n";
      if(shift.getStart() >= timeWindows[tw].first && end <= timeWindows[tw].second
              /*&& shift.getStart() <= timeWindows[tw].second && end >= timeWindows[tw].first*/
        && shift.getStart()<= end && timeWindows[tw].first<=timeWindows[tw].second){
            flag = false;
      }

    }
    
    if(flag){
       COUT<<"DRI08 NOT FEASIBLE!!!!";
       int a; CIN>>a;
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
//	  COUT<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<end1<<" \n";

      if(  not(  shift2.getStart() >= end1 //or shift1.getStart()  >= end2
                 )){
         COUT<<"TL01 NOT FEASIBLE!!!!";
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
    //	  COUT<<"CHECK: "<<solution.getShifts().at(s1).getIndex()<< " " << solution.getShifts().at(s2).getIndex()<<" "<<end1<<" \n";

          if(  not(  shift2.getStart() >= end1 //or shift1.getStart()  >= end2
                     )){
             COUT<<"TL01 NOT FEASIBLE!!!!";
             int a; CIN>>a;
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
    
//    COUT<<solution.getShifts()[s].getTrailer().getIndex()<<" "<<solution.getShifts()[s].getDriver().getIndex()<<" "
//    <<this->drivers[solution.getShifts()[s].getDriver().getIndex()].getTrailer()<<"\n";
    if( not(solution.getShifts()[s].getTrailer()==this->drivers[solution.getShifts()[s].getDriver()].getTrailer() ) ){
       COUT<<"TL03 NOT FEASIBLE!!!!";
       int a; CIN>>a;
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

  double unfeasibilityCounter = 1/EPSILON;

  vector<vector<double> > tankQuantities(this->customers.size());
  for(int p=0; p<this->customers.size();p++){
    if(p>1){
      vector<double> horizon(this->customers[p].getForecast().size()/**60*/,0.0);
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
             COUT<<"SHI11 NOT FEASIBLE!!!!";
                COUT<<"\n"<<s<<" "<<operations[o].getPoint()<<" "<<operations[o].getQuantity()<<" ";
                int a; CIN>>a;
//                unfeasibilityCounter = operations[o].getArrival();
              /*exit(0);*/
//              return true;
          }
          tankQuantities[operations[o].getPoint()][(int)(operations[o].getArrival()/60)] += operations[o].getQuantity();
        }
        else if(operations[o].getPoint() == 1)
          if(not (operations[o].getQuantity() <= EPSILON)){
             COUT<<"SHI11 NOT FEASIBLE!!!!";
             int a; CIN>>a;
//             unfeasibilityCounter = operations[o].getArrival();
              /*exit(0);*/
//              return true;
          }
    }
  }

  bool breakFlag = false;
  for(int p=0; p<this->customers.size(); p++){
      vector<double> fore = this->customers[p].getForecast();
    if(p>1){
      for(int f=0; f<fore.size()/**60*/; f++){
 //       if(f%60 == 0)
        tankQuantities[p][f] -= fore[(int)f/*/60*/];
        if(f>=1)
          tankQuantities[p][f] += tankQuantities[p][f-1];

        if( not(  (this->customers[p].getCapacity() - tankQuantities[p][f] >= -EPSILON or abs(this->customers[p].getCapacity() - tankQuantities[p][f]) <= EPSILON)
          && tankQuantities[p][f] >= -EPSILON
          && ( tankQuantities[p][f] >= this->customers[p].getSafetyLevel() or abs(tankQuantities[p][f] - this->customers[p].getSafetyLevel()) <= EPSILON)
                   )
             ){
//               COUT<<"DYN01 QS02 NOT FEASIBLE!!!!";
                  if(abs(f) < unfeasibilityCounter)
                    unfeasibilityCounter = abs(f);

                  breakFlag = true;
                  break;
            }

      }
      /* if(breakFlag){
           break;
       }*/
    }
  }
  if(unfeasibilityCounter >= 1/EPSILON - EPSILON)
      return 0;
  else
  return unfeasibilityCounter;
}


double Instance::dyn01partial(irpSolution solution, bool feasibility, unsigned int partialHorizon){


  partialHorizon = (unsigned int)(partialHorizon/60) + 1;
  COUT<<"PARTIAL HORIZON: "<<partialHorizon<<"\n";
  double unfeasibilityCounter = 1/EPSILON;

  vector<vector<double> > tankQuantities(this->customers.size());
  for(int p=0; p<this->customers.size();p++){
    if(p>1){
      vector<double> horizon(this->customers[p].getForecast().size()/**60*/,0.0);
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
             COUT<<"SHI11 NOT FEASIBLE!!!!";
                COUT<<"\n"<<s<<" "<<operations[o].getPoint()<<" "<<operations[o].getQuantity()<<" ";
                int a; CIN>>a;
//                unfeasibilityCounter = operations[o].getArrival();
              /*exit(0);*/
//              return true;
          }
          tankQuantities[operations[o].getPoint()][(int)(operations[o].getArrival()/60)] += operations[o].getQuantity();
        }
        else if(operations[o].getPoint() == 1)
          if(not (operations[o].getQuantity() <= EPSILON)){
             COUT<<"SHI11 NOT FEASIBLE!!!!";
             int a; CIN>>a;
//             unfeasibilityCounter = operations[o].getArrival();
              /*exit(0);*/
//              return true;
          }
    }
  }

  for(int p=0; p<this->customers.size(); p++){
    if(p>1){
      for(int f=0; f<partialHorizon/**60*/; f++){
 //       if(f%60 == 0)
        tankQuantities[p][f] -= this->customers[p].getForecast()[(int)f/*/60*/];
        if(f>=1)
          tankQuantities[p][f] += tankQuantities[p][f-1];

        if( not(  (this->customers[p].getCapacity() - tankQuantities[p][f] >= -EPSILON or abs(this->customers[p].getCapacity() - tankQuantities[p][f]) <= EPSILON)
          && tankQuantities[p][f] >= -EPSILON
          && ( tankQuantities[p][f] >= this->customers[p].getSafetyLevel() or abs(tankQuantities[p][f] - this->customers[p].getSafetyLevel()) <= EPSILON)
                   )
             ){
//               COUT<<"DYN01 QS02 NOT FEASIBLE!!!!";
                  if(abs(f) < unfeasibilityCounter)
                    unfeasibilityCounter = abs(f);
            }
      }
    }
  }
  if(unfeasibilityCounter >= 1/EPSILON - EPSILON)
      return 0;
  else
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
        COUT<<"SHI03 NOT FEASIBLE!!!!";
        int a; CIN>>a;
         /*exit(0);*/
         return true;
        }

        for(int o=1; o<operations.size(); o++){
          departure = operations[o-1].getArrival() + this->customers[operations[o-1].getPoint()].getSetupTime();
    //      if(s==1)
    //      COUT<<operations[o].getArrival()<<" "<<departure + this->timeMatrices[operations[o-1].getPoint()][operations[o].getPoint()]<<"\n";
          if(not(  operations[o].getArrival() >= departure + this->timeMatrices[operations[o-1].getPoint()][operations[o].getPoint()] )){
         COUT<<"SHI02 NOT FEASIBLE!!!!\n";
         int a; CIN>>a;
         COUT<<s<<" "<<o<<" "<<departure + this->timeMatrices[operations[o-1].getPoint()][operations[o].getPoint()]<<" "<<operations[o].getArrival()<<"\n";
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

      int a; if(operations[o].getPoint() == 0)CIN>>a;

      //Every trailer is allowed at the source location
      if(operations[o].getPoint()==1)
        flag=true;
      for(int t=0; t<allowedTrailers.size(); t++ ){
	
        if( shift.getTrailer()==allowedTrailers[t] )
          flag = true;
      }
      if(not(flag)){
         COUT<<"SHI05 NOT FEASIBLE!!!!\n";
         COUT<<s<<" "<<o<<"\n";
         CIN>>a;
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
            COUT<<"SHI06 NOT FEASIBLE!!!!  ";
             COUT<<s<<" "<<o<<"     "<<oper<<"    "<<trailerQuantities[shift.getTrailer()][oper]<<"\n";
             int a; CIN>>a;
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

    //Time to check all constraints is no relevant
    double feas = this->dyn01(solution, false);
    if(not feas){
        if(
 /*             this->dri01(solution) ||
              this->dri03(solution) ||
              this->dri08(solution) ||
              this->tl01(solution)  ||
              this->tl03(solution)  ||*/
              feas /*||
              this->shi02(solution) ||
              this->shi05(solution) ||
              this->shi06(solution)*/
                      )
            return feas;
        else
            return 0;
    }
    else
        return feas;

}

double Instance::computeObjective(irpSolution &solution){
  
  double totalQuantity = 0.0;
  double shiftQuantity, shiftDistance, shiftTime, shiftCost;
  shiftQuantity = shiftDistance = shiftTime = shiftCost = 0.0;
  double LR = 0.0;
  for(int s=0; s<solution.getShifts().size(); s++){

    double travelDist = 0.0;
    /*shiftQuantity = 0.0;*/
    double cost = 0.0;

    Shift shift = solution.getShifts()[s];
    vector<Operation> operations = shift.getOperations();

    travelDist += this->distMatrices[0][operations.front().getPoint()];

    for(int o=0; o<operations.size(); o++){
      
      if(operations[o].getQuantity()>0){
        totalQuantity += operations[o].getQuantity();
        shiftQuantity += operations[o].getQuantity();
      }
      
      if(o>0)
        travelDist += this->distMatrices[operations[o-1].getPoint()][operations[o].getPoint()];
    }
    travelDist+=this->distMatrices[operations.back().getPoint()][0];


    double end = operations.back().getArrival() + this->customers[operations.back().getPoint()].getSetupTime() + this->timeMatrices[operations.back().getPoint()][0];
    cost = this->trailers[shift.getTrailer()].getDistanceCost()*travelDist + this->drivers[shift.getDriver()].getTimeCost()*(end - shift.getStart());
    LR+=cost;

    shiftDistance += travelDist;
    shiftTime += (end - shift.getStart());
    shiftCost += cost;

    shift.setShiftQuantity(shiftQuantity);
    shift.setShiftDistance(shiftDistance);
    shift.setShiftTime(shiftTime);
    shift.setShiftCost(shiftCost);
    solution.setShift(shift, s);
  }


  solution.setTotalQuantity(totalQuantity);
  solution.setCost(LR);
  LR/=totalQuantity;
//  COUT<<"COMPUTE: "<<solution.getShifts()[2].getShiftQuantity()<<" "<<solution.getShifts()[2].getShiftDistance()<<" "<<solution.getShifts()[2].getShiftTime()<<"\n";
  return LR;
}


double Instance::computeDeltaObjective(irpSolution oldSolution, irpSolution newSolution, unsigned int shiftIndex){

    Shift oldShift = oldSolution.getShifts()[shiftIndex];
    double oldQuantity = oldSolution.getShifts()[shiftIndex].getShiftQuantity();
    double oldDistance = oldSolution.getShifts()[shiftIndex].getShiftDistance();
    double oldTime = (oldShift.getOperations().back().getArrival()
                      + this->customers[oldShift.getOperations().back().getPoint()].getSetupTime()
                      + this->timeMatrices[oldShift.getOperations().back().getPoint()][0]
                      - oldShift.getStart());
    double oldShiftObjective = oldDistance * this->trailers[oldShift.getTrailer()].getDistanceCost() +
                               oldTime * this->drivers[oldShift.getDriver()].getTimeCost();


    Shift newShift = newSolution.getShifts()[shiftIndex];
    double newQuantity = 0.0;
    double newDistance = 0.0;
    double newTime = (newShift.getOperations().back().getArrival()
                      + this->customers[newShift.getOperations().back().getPoint()].getSetupTime()
                      + this->timeMatrices[newShift.getOperations().back().getPoint()][0]
                      - newShift.getStart());

    vector<Operation> newOperations = newShift.getOperations();
    newDistance += this->distMatrices[0][newOperations.front().getPoint()];
    for(int o=0; o<newOperations.size(); o++){
        if(newOperations[o].getQuantity() > EPSILON)
            newQuantity += newOperations[o].getQuantity();
        if(o > 0)
            newDistance+= this->distMatrices[newOperations[o-1].getPoint()][newOperations[o].getPoint()];
    }
    newDistance += this->distMatrices[newOperations.back().getPoint()][0];

    double newShiftObjective = newDistance * this->trailers[newShift.getTrailer()].getDistanceCost() +
                               newTime * this->drivers[newShift.getDriver()].getTimeCost();


    newShiftObjective /= oldSolution.getTotalQuantity();
    oldShiftObjective /= (oldSolution.getTotalQuantity() - oldQuantity + newQuantity);


    double newObjectiveValue = oldSolution.getCost()/oldSolution.getTotalQuantity();
    newObjectiveValue -= oldShiftObjective;
    newObjectiveValue *= oldSolution.getTotalQuantity();
    newObjectiveValue /= (oldSolution.getTotalQuantity() - oldQuantity + newQuantity);
    newObjectiveValue += newShiftObjective;

//    COUT<<"\nINSIDE DELTA: "<<oldShiftObjective<<" "<<newShiftObjective<<" "<<oldQuantity<<" "<<newQuantity<<"\n";

    return newObjectiveValue;

}

double Instance::computePartialObjective(irpSolution oldSolution, irpSolution newSolution, unsigned int shiftIndex){

    double oldQuantity, newQuantity;
    double oldDistance, newDistance;
    double oldTime, newTime;
    double oldShiftObjective;

    oldQuantity = newQuantity = oldDistance = newDistance = oldTime = newTime = oldShiftObjective = 0.0;

    if(shiftIndex > 0){
        oldQuantity += oldSolution.getShifts()[shiftIndex-1].getShiftQuantity();
//         COUT<<"OLDQ: "<<oldSolution.getShifts()[shiftIndex-1].getShiftQuantity()<<"\n";
        oldDistance = oldSolution.getShifts()[shiftIndex-1].getShiftDistance();
        oldTime = oldSolution.getShifts()[shiftIndex-1].getShiftTime();
        oldShiftObjective = oldSolution.getShifts()[shiftIndex-1].getShiftCost();
    }


//    COUT<<"OLD QUANTITY: "<<oldQuantity<<" "<<oldDistance<<" "<<oldTime<<" "<<oldShiftObjective<<"\n";

    double newShiftObjective = 0.0;
    for(int s=shiftIndex; s<newSolution.getShifts().size(); s++){
        newDistance = newTime = 0.0;
        Shift newShift = newSolution.getShifts()[s];
        vector<Operation> newOperations = newShift.getOperations();
        newDistance += this->distMatrices[0][newOperations.front().getPoint()];
        for(int o=0; o<newOperations.size(); o++){
            if(o > 0)
                newDistance += this->distMatrices[newOperations[o-1].getPoint()][newOperations[o].getPoint()];
            if(newOperations[o].getPoint() != 1)
                newQuantity += newOperations[o].getQuantity();
        }
        newDistance += this->distMatrices[newOperations.back().getPoint()][0];
        newTime = (newOperations.back().getArrival() + this->customers[newOperations.back().getPoint()].getSetupTime() +  this->timeMatrices[newOperations.back().getPoint()][0] - newShift.getStart());


        newShiftObjective += (newDistance * this->trailers[newShift.getTrailer()].getDistanceCost() +
                                   newTime * this->drivers[newShift.getDriver()].getTimeCost());
    }

    double actualQuantity = oldQuantity + newQuantity;
    double newObjective = (oldShiftObjective + newShiftObjective)/actualQuantity;

    return newObjective;
}

double Instance::computePartialDeltaObjective(irpSolution oldSolution, irpSolution newSolution, unsigned int shiftIndex){


    double p = this->computePartialObjective(oldSolution, newSolution, shiftIndex);

//    clock_t begin = clock();
    double feas = this->checkFeasibility(newSolution, false);
    /*
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    COUT<<"INNER ELAPSED: "<<elapsed_secs<<"\n";
    */

    if(not(feas < EPSILON and feas > -EPSILON)){
        double factor = feas/(this->horizon*60);
        p /=  (factor * FEASIBILITY_PENALTY);
        p += 1.0;
        COUT<<"\n NOT FEASIBLE! \n";
    }

    return p;
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

 /* vector<double> addP(priorities.size(), 0.0);
  for(int p=0; p<addP.size(); p++)
      addP[p] = 100 * EPSILON * this->drivers[priorities[p].second].getTimeCost()
              + 100 * EPSILON * this->trailers[this->drivers[priorities[p].second].getTrailer()].getDistanceCost()
              - 100 * EPSILON * this->trailers[this->drivers[priorities[p].second].getTrailer()].getInitialQuantity()/this->maxInitialQuantity
              - 100 * EPSILON * this->trailers[this->drivers[priorities[p].second].getTrailer()].getCapacity()/this->maxCapacity;

 /* COUT<<setprecision(15);
  for(int p=0; p<addP.size(); p++)
      COUT<<(double)((double)priorities[p].first.first + addP[p])<<"    ";
*/
  pair< pair<unsigned int,unsigned int>, unsigned int > app;
  double a;
    for (int i=0; i<priorities.size(); i++){
        for (int j=i+1; j<priorities.size(); j++){

            if ((double)((double)priorities[i].first.first /*+ addP[i]*/) > (double)((double)priorities[j].first.first /*+ addP[j]*/)){
                app =  priorities[i];
                priorities[i] = priorities[j];
                priorities[j] = app;

               /* a =  addP[i];
                addP[i] = addP[j];
                addP[j] = a;*/
            }
        }
    }
}



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
                                                 vector<Operation> &operations,
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
        if(customerList[succ] != 1){
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
     //max driving duration constraint
     cumTime + this->timeMatrices[customerList[prec]][customerList[succ]] + this->customers[customerList[succ]].getSetupTime() + this->timeMatrices[customerList[succ]][0] <= this->drivers[timeWindows[0].second].getMaxDrivingDuration()
     //shift within time window
     and timeWindows[0].first.first + cumTime + this->timeMatrices[customerList[prec]][customerList[succ]] + this->customers[customerList[succ]].getSetupTime() + this->timeMatrices[customerList[succ]][0] <= timeWindows[0].first.second
     //shifts with same driver should be distanciated by a minimum intershift duration
     and (timeWindows[twD].first.first - (timeWindows[0].first.first + cumTime + this->timeMatrices[customerList[prec]][customerList[succ]] + this->customers[customerList[succ]].getSetupTime() + this->timeMatrices[customerList[succ]][0]) >= drivers[timeWindows[0].second].getMinInterSHIFTDURATION()
       or lastDShift   )
     //shift using the same trailer cannot overlap
     and (timeWindows[twT].first.second <= timeWindows[0].first.first or timeWindows[twT].first.first >= timeWindows[0].first.second
      or lastTShift    )
     //customer site should be accessible by the shift's trailer
     and (allowedTrailerFlag or customerList[succ]==1)
     )


     or (tankQuantities[customerList[succ]] >= this->customers[customerList[succ]].getCapacity()-EPSILON and customerList[succ] != 1)
     or (lastOperations[customerList[succ]] >= time + this->timeMatrices[customerList[prec]][customerList[succ]]  and customerList[succ] != 1)
     //demand for the customer is completely satisfied
     or (horizon >= -EPSILON and customerList[succ] != 1)
  ){

     //it is not possible to go to the source to refuel, shift ends
     if(customerList[succ] == 1)
         return operations;
     customerList.erase(customerList.begin() + succ);

     //all customers have been served or cannot be served, shift ends
     if(customerList.size() <= 1 )
       return operations;


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

    // time at which last operation for this customer has been performed
    if(lastOperations[customer] < time + this->timeMatrices[customerList[prec]][customer])
        lastOperations[customer] = time + this->timeMatrices[customerList[prec]][customer];

    operation.setPoint(customer);


//    COUT<<"         Tank: "<<tankQuantities[customer]<<"   trail:"<<trailerQuantities[trailer]<<"   horizon: "<<horizon<<" "<<this->customers[customer].getCapacity()<< "\n";


    if(customer == 1){
      operation.setQuantity(-(this->trailers[trailer].getCapacity() - trailerQuantities[trailer]));
      pastQuantities[customer] += operation.getQuantity();
      trailerQuantities[trailer] -= /*-*/ operation.getQuantity();
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
      vector<double> forecast = this->customers[c].getForecast();
      double totForecast = 0;
      bool flag = true;

      double deliveredQuantity = deliveredQuantities[c];

      unsigned int currentPoint = customer;

      double totalDistance = 0.0;
      for(int d=0; d<this->customers.size(); d++)
          totalDistance += this->timeMatrices[currentPoint][d];
      totalDistance = *max_element(this->timeMatrices[currentPoint].begin(), this->timeMatrices[currentPoint].end());
      double totalForecast = totalForecasts[c];

      unsigned int f = 0;
      double timeFactor, quantityFactor, distanceFactor;
      while(f < forecast.size() and flag){
        totForecast += forecast[f];
        if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
          timeFactor = ((double)f/forecast.size()) * timeWeight/2;
          quantityFactor = (abs(totalForecast - deliveredQuantity)/totalForecast) * quantityWeight/2;
          distanceFactor = 100.0 * EPSILON * ((double)this->timeMatrices[currentPoint][c]/totalDistance)*ties;
          urgency[c] = 1.0 + timeFactor + quantityFactor + distanceFactor;
          flag = false;

  //        COUT<<"Fs: "<<(double)this->timeMatrices[currentPoint][c]<<" "<<totalDistance<<"\n";
  //        COUT<<"Fs: "<<c<<" "<<timeFactor<<" "<<quantityFactor<<" "<<distanceFactor<<" "<<urgency[c]<<"\n\n";
//            int a;CIN>>a;
//          COUT<<c<<" "<<distanceFactor<<" "<<urgency[c]<<"\n";
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

//    COUT<<"     Operation: "<<operation.getPoint()<<"   arr: "<<operation.getArrival()<<"   q: "<<operation.getQuantity()<<"\n";
    /*
    for(int c=0; c<customerList.size(); c++)
      COUT<<customerList[c]<<" ";
    COUT<<"\n";

    for(int u=0; u<urgency.size(); u++)
      COUT<<urgency[u]<<" ";
    COUT<<"\n";
    */

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
        /////ADDED
        while(solution.getShifts().back().getOperations().size() == 1 and solution.getShifts().back().getOperations().front().getPoint() == 1){
        /////
            vector<Shift> shifts = solution.getShifts();
            shifts.erase(shifts.end());
            solution.setShifts(shifts);
       }

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
        totalDistance = *max_element(this->timeMatrices[currentPoint].begin(), this->timeMatrices[currentPoint].end());
        double totalForecast = totalForecasts[c];

        unsigned int f = 0;
        double timeFactor, quantityFactor, distanceFactor;
        while(f < forecast.size() and flag){
          totForecast += forecast[f];
          if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
              timeFactor = ((double)f/forecast.size()) * timeWeight/2;
              quantityFactor = (abs(totalForecast - deliveredQuantity)/totalForecast) * quantityWeight/2;
              distanceFactor = 100.0 * EPSILON * ((double)this->timeMatrices[currentPoint][c]/(totalDistance))*ties;
              urgency[c] = 1.0 + timeFactor + quantityFactor + distanceFactor;
              flag = false;

      //        COUT<<"Fs: "<<(double)this->timeMatrices[currentPoint][c]<<" "<<totalDistance<<"\n";
      //        COUT<<"Fs: "<<c<<" "<<timeFactor<<" "<<quantityFactor<<" "<<distanceFactor<<" "<<urgency[c]<<"\n\n";
    //            int a;CIN>>a;
    //          COUT<<c<<" "<<totalForecast - this->customers[c].getInitialTankQuantity() + this->customers[c].getSafetyLevel()<<" "<<urgency[c]<<"\n";
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
      /////ADDED
      customerList.insert(customerList.begin(), 1);
      /////
      customerList.insert(customerList.begin(), 0);

/*
      for(int c=0; c<customerList.size(); c++)
        COUT<<customerList[c]<<";
      COUT<<"\n";

      for(int u=0; u<urgency.size(); u++)
        COUT<<urgency[u]<<" ";
      COUT<<"\n";
*/
//      int a;CIN>>a;
//      COUT<<"SHIFT: "<<shift.getIndex()<<" start: "<<shift.getStart()<<"  \n";
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
         /////ADDED
         if(not(operations.size()==1 and operations.front().getPoint() == 1 and abs(operations.front().getQuantity()) <= EPSILON)){
         /////
             /////ADDED
             //Refuelling at the end of a shift, when at the begginning of a shift a refuel is performed, is useless
             if(operations.size()>1 and operations.back().getPoint() == 1){
                 trailerQuantities[shift.getTrailer()] += operations.back().getQuantity();
                 tankQuantities[operations.back().getPoint()] -= operations.back().getQuantity();
                 operations.erase(operations.end());
             }
             /////

             shift.setOperations(operations);
             vector <Shift> shifts = solution.getShifts();
             shifts.push_back(shift);
             solution.setShifts(shifts);
         }
      }


      timeWindows.erase(timeWindows.begin());
      return recursiveRandomSolution(solution, time, horizonQuantities, tankQuantities, trailerQuantities,
                                     pastQuantities, lastOperations, timeWindows,
                                     timeWeight, quantityWeight, ties,
                                     deliveredQuantities, totalDistance, totalForecasts);
    }
}

irpSolution Instance::backTrackingRandomSolution(double timeWeight, double quantityWeight, double ties){

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
      vector<double> horizon(this->horizon,0.0);
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
      totalForecasts[c] = abs(horizonQuantities[c].back());
     // for(int f=0; f<this->customers[c].getForecast().size(); f++)
     //     totalForecasts[c] += this->customers[c].getForecast()[f];


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
                                      double refuelRatio, double deliveredQuantityRatio, bool originalFlag){

   deliveredQuantityRatio = 1.0 - deliveredQuantityRatio;

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

       COUT<<"SHIFT: "<<shift.getIndex()<<" start: "<<shift.getStart()<<" "<<shift.getDriver()<<" "<<shift.getTrailer()<<"\n";
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
    /*
       /////MODIFIED
           prec = 0;
           succ = 1;
       /////MODIFIED
       /// */
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

//           COUT<<"         Tank: "<<this->actualQuantity[succ][operation.getArrival()]<<"   trail:"<<trailerQuantities[trailer]<<" "<<this->customers[succ].getCapacity()<<"       \n";
           double quantity = 0.0;

           if(refuelFlag[trailer] == true or succ == 1){
               quantity = this->trailers[trailer].getCapacity() - trailerQuantities[trailer];
               trailerQuantities[trailer] += quantity;
               refuelFlag[trailer] = false;

               operation.setQuantity(-quantity);

               COUT<<"     Operation: "<<operation.getPoint()<<"   arr: "<<operation.getArrival()<<"   q: "<<operation.getQuantity()<<"\n";
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
               if(trailerQuantities[trailer] <= fabs(this->horizons[succ].back())
                       and this->horizons[succ].back() < 0){
                   double residualQuantity = this->actualQuantity[succ][(int)(operation.getArrival())];
                   if(residualQuantity < 0)
                       residualQuantity = 0;
                   else if(residualQuantity > this->customers[succ].getCapacity())
                       residualQuantity = this->customers[succ].getCapacity();


                       if(trailerQuantities[trailer] <= this->customers[succ].getCapacity() - residualQuantity){
                           quantity = trailerQuantities[trailer];

                           /////check if the quantity would make the tank exceed due to future operations, if so break
                           for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                               if(this->actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                                  quantity -= (this->actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                                  horizonFlag = true;
                               }
                           }
                           if(horizonFlag)
                               break;
                           //////

                           /////Delivered quantity factor
                           double factor = 1.0;// - residualQuantity/this->customers[succ].getCapacity()
                           if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                            quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                           /////

                           trailerQuantities[trailer] -= quantity;
                           tankQuantities[succ] += quantity;

                       }
                       else{
                           quantity = (this->customers[succ].getCapacity() - residualQuantity);

                           /////check if the quantity would make the tank exceed due to future operations, if so break
                           for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                               if(this->actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                                  quantity -= (this->actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                                  horizonFlag = true;
                               }
                           }
                           if(horizonFlag)
                               break;
                           /////

                           /////Delivered quantity factor
                           double factor = 1.0;// - residualQuantity/this->customers[succ].getCapacity();
                           if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                               quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                           /////

                           trailerQuantities[trailer] -= quantity;
                           tankQuantities[succ] += quantity;
                       }
               }
           else{
               double residualQuantity = this->actualQuantity[succ][(int)(operation.getArrival())];
               if(residualQuantity < 0)
                   residualQuantity = 0;
               else if(residualQuantity > this->customers[succ].getCapacity())
                   residualQuantity = this->customers[succ].getCapacity();

               if(fabs(this->horizons[succ].back()) <= this->customers[succ].getCapacity() - residualQuantity
                       and this->horizons[succ].back() < 0){
                   quantity = fabs(this->horizons[succ].back());

                   /////check if the quantity would make the tank exceed due to future operations, if so break
                   for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                       if(this->actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                          quantity -= (this->actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                          horizonFlag = true;
                       }
                   }
                   if(horizonFlag)
                       break;
                   /////

                   /////Delivered quantity factor
                   double factor = 1.0;
                   if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                      quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                   /////

                   trailerQuantities[trailer] -= quantity;
                   tankQuantities[succ] += quantity;
               }
               else if(fabs(this->horizons[succ].back()) > this->customers[succ].getCapacity() - residualQuantity
                       and this->horizons[succ].back() < 0){
                   quantity = (this->customers[succ].getCapacity() - residualQuantity);

                   /////check if the quantity would make the tank exceed due to future operations, if so break
                   for(int ff=operation.getArrival(); ff < this->horizon*60; ff++){
                       if(this->actualQuantity[succ][ff] + quantity > this->customers[succ].getCapacity()){
                          quantity -= (this->actualQuantity[succ][ff] + quantity - this->customers[succ].getCapacity());
                          horizonFlag = true;
                        }
                   }
                   if(horizonFlag)
                       break;
                   /////

                   /////Delivered quantity factor
                   double factor = 1.0;// - residualQuantity/this->customers[succ].getCapacity();
                   if(not(quantity <= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio))
                        quantity -= this->trailers[trailer].getCapacity() * factor * deliveredQuantityRatio;
                   /////

                   trailerQuantities[trailer] -= quantity;
                   tankQuantities[succ] += quantity;

               }
               else
                   quantity = 0.0;
           }

           if(horizonFlag or (quantity <= EPSILON and operation.getPoint()!=1 /*and originalFlag*/))
                break;
           if(trailerQuantities[trailer] <= EPSILON + this->trailers[trailer].getCapacity() * refuelRatio){
             refuelFlag[trailer] = true;
           }

           for(int ff = (int)(operation.getArrival()); ff < this->customers[succ].getForecast().size()*60; ff++){
               this->horizons[succ][ff] += quantity;
               this->actualQuantity[succ][ff] += quantity;
           }

           operation.setQuantity(quantity);

           COUT<<"     Operation: "<<operation.getPoint()<<"   arr: "<<operation.getArrival()<<"   q: "<<operation.getQuantity()<<"\n";
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
       timeWindows.erase(timeWindows.begin());
       shift.setOperations(operations);
       if(operations.size() != 0)
        shifts.push_back(shift);
       else
           index--;

       index++;
   }


   vector<Shift> newShifts;
   unsigned int size = shifts.size();
   for(int s=size-1; s>=0; s--){
       Shift newShift = shifts[s];
       vector<Operation> newOperations = shifts[s].getOperations();
       unsigned int osize = shifts[s].getOperations().size();
       for(int o=osize-1; o>=0; o--)
           if(shifts[s].getOperations()[o].getQuantity() < EPSILON && shifts[s].getOperations()[o].getQuantity() > -EPSILON )
               newOperations.erase(newOperations.end());
           else
               break;

       newShift.setOperations(newOperations);
       if(newShift.getOperations().size() > 0)
        shifts[s] = newShift;
       else
           shifts.erase(shifts.end());
   }
//   shifts = newShifts;

      solution.setShifts(shifts);

   for(int s=0; s<solution.getShifts().size(); s++){
       COUT<<"SHIFT R: "<<solution.getShifts()[s].getIndex()<<" "<<solution.getShifts()[s].getStart()<<" "<<solution.getShifts()[s].getDriver()<<"\n";
       for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++)
           COUT<<"     "<<solution.getShifts()[s].getOperations()[o].getPoint()<<" "
                 <<solution.getShifts()[s].getOperations()[o].getArrival()<<" "
                   <<solution.getShifts()[s].getOperations()[o].getQuantity()<<"\n";
   }
   return solution;
}

unsigned int Instance::computeTime(unsigned int succ, unsigned int nextSucc){

    vector<unsigned int> travellingTimes = this->timeMatrices[succ];

//    COUT<<"\n"<<succ<<"times: \n";
    for(int c=0; c<travellingTimes.size(); c++){
//        COUT<<travellingTimes[c]<<" ";
        travellingTimes[c] += this->customers[c].getSetupTime();
//        COUT<<travellingTimes[c]<<" ";
        travellingTimes[c] += this->timeMatrices[c][0];
//        COUT<<travellingTimes[c]<<"\n";
    }

//    COUT<<"MIN: "<<*std::min_element(travellingTimes.begin(), travellingTimes.end())<<"\n";
    return *std::min_element(travellingTimes.begin(), travellingTimes.end());

}

void Instance::computeUrgency2(vector<unsigned int> &customerList, vector<double> &urgency,
                               Shift shift, vector<double> deliveredQuantities,
                               double costWeight, double demandWeight,
                               unsigned int prec){

    vector <vector <double> > cc = this->contributeMatrixesValues.find(shift.getDriver())->second;
    vector<double> v = cc[prec];
    for(int i=0; i<v.size(); i++)
        if(v[i] >= DBL_MAX- EPSILON)
            v[i] = 0.0;

    double maxC  = *max_element(v.begin()+2, v.end());
    double demandFactor = 0.0;
    double costFactor;
    for(int c=0; c<urgency.size(); c++){
        if(c > 1){
            costFactor = ((double)v[c]/maxC);

            COUT<<"COST: "<<costFactor<<" "<<v[c]<<" "<<maxC<<"\n";
            Customer customer = this->customers[c];
            vector<double> forecast = customer.getForecast();
            double deliveredQuantity = deliveredQuantities[c];
            double totForecast = 0;
            unsigned int f = 0;
            bool flag = true;
            while(f < forecast.size() and flag){
              totForecast += forecast[f];
              if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
                  demandFactor = (double)f/forecast.size();
                  flag = false;
                  COUT<<"DEMAND: "<<demandFactor<<" "<<f<<"\n";
              }
              f++;
            }


            if(not flag)
              urgency[c] = (costFactor*costWeight + demandFactor*demandWeight)/2;
            else
                urgency[c] = 1.0 + costFactor*costWeight;
        }
        else
            urgency[c] = DBL_MAX;

    }

}

irpSolution Instance::constructSolution2(irpSolution initialSolution, double costWeight, double demandWeight, int ties,
                                        double servingRatio, double refuelRatio,
                                        unsigned int noCustomer,
                                        unsigned int maxShifts){

    irpSolution solution = initialSolution;

    vector<Shift> shifts = solution.getShifts();

    if(shifts.size() > maxShifts)
        return solution;

    vector<double> trailerQuantities;
    for(int t=0; t<this->trailers.size(); t++)
        trailerQuantities.push_back(this->trailers[t].getInitialQuantity());

    vector<double> tankQuantities;
    for(int c=0; c<this->customers.size(); c++)
        tankQuantities.push_back(this->customers[c].getInitialTankQuantity());


    vector<double> deliveredQuantities(this->customers.size(), 0.0);

    vector<Operation> lastOperations(this->customers.size(), Operation());
    vector<unsigned int> lastOperationsTime(this->customers.size(), 0);
    unsigned int lastTime = 0;
    for(int s=0; s<shifts.size(); s++){
        Shift shift = shifts[s];
        for(int o=0; o<shift.getOperations().size(); o++){
            Operation operation = shift.getOperations()[o];
            trailerQuantities[shift.getTrailer()] -= operation.getQuantity();
            tankQuantities[operation.getPoint()] += operation.getQuantity();
            deliveredQuantities[operation.getPoint()] += operation.getQuantity();

            if(operation.getArrival() > lastOperationsTime[operation.getPoint()]){
                lastOperationsTime[operation.getPoint()] = operation.getArrival();
                lastOperations[operation.getPoint()] = operation;
            }

        }
        if(shift.getStart()  > lastTime)
            lastTime = shift.getStart();
    }

    vector< pair< pair<unsigned int,unsigned int>, unsigned int > > lastTimeWindows;

    for(int t=0; t<this->timeWindows.size(); t++){
      if(this->timeWindows[t].first.first >= lastTime){

          pair< pair<unsigned int,unsigned int>, unsigned int > tw;
          tw.first = this->timeWindows[t].first;
          tw.second = this->timeWindows[t].second;
          lastTimeWindows.push_back(tw);

          int numberOfDrivers = this->drivers.size();

          if(numberOfDrivers > shifts.size())
              numberOfDrivers -= (numberOfDrivers - shifts.size());

          for(int d=0; d<numberOfDrivers; d++){
              Operation lastOperation = shifts[/*s*/shifts.size()-1-d].getOperations().back();
              int lastOperationTime = lastOperation.getArrival()
                      + this->customers[lastOperation.getPoint()].getSetupTime()
                      + this->timeMatrices[lastOperation.getPoint()][0];
              if( (abs(lastOperationTime - (int)tw.first.first) < this->drivers[tw.second].getMinInterSHIFTDURATION() or shifts[/*s*/initialSolution.getShifts().size()-1-d].getStart()==tw.first.first)
                     and (int)tw.second == (int)shifts[/*s*/shifts.size()-1-d].getDriver() )
                  lastTimeWindows.erase(lastTimeWindows.end());
          }
      }
    }

    if(lastTimeWindows.size() == 0)
        return solution;

    unsigned int time, index = 0;
    if(shifts.size() != 0){
        time = lastTimeWindows.front().first.first;
        index = shifts.back().getIndex() + 1;
    }


    vector<bool> refuelFlags(this->trailers.size(), false);
    for(int t=0; t<this->trailers.size(); t++)
        if(trailerQuantities[t] < EPSILON)
            refuelFlags[t] = true;
//    vector<Shift> shifts = solution.getShifts();

    while(lastTimeWindows.size() > 0)  {

        Shift shift;
        shift.setIndex(index);
        shift.setDriver(lastTimeWindows.front().second);
        shift.setTrailer(this->drivers[shift.getDriver()].getTrailer());
        shift.setStart(lastTimeWindows.front().first.first);

//        COUT<<"SHIFT C: "<<shift.getIndex()<<" "<<shift.getStart()<<" "<<shift.getDriver()<<" "<<shift.getTrailer()<<"\n";

        time = shift.getStart();

        unsigned int prec;
        unsigned int succ;

        prec = 0;
        succ = 1;

        unsigned int cumTime = 0;

        unsigned int twD=0;
        bool lastDShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(lastTimeWindows[t].second == shift.getDriver()){
            twD = t;
            lastDShift = false;
            break;
          }
        unsigned int twT=0;
        bool lastTShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(this->drivers[lastTimeWindows[t].second].getTrailer() == this->drivers[shift.getDriver()].getTrailer()){
            twT = t;
            lastTShift = false;
            break;
          }
          bool allowedTrailerFlag = false;
          for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
            if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[shift.getDriver()].getTrailer())
                allowedTrailerFlag = true;


          vector<unsigned int> customerList;
          for(int c=0; c<this->customers.size(); c++)
              customerList.push_back(c);

          vector<double> urgency(this->customers.size(), DBL_MAX);
////////////////////////////////////

          this->computeUrgency2(customerList, urgency,
                                         shift, deliveredQuantities,
                                         costWeight, demandWeight, 0);

          ////////////////////////////////



          for(int u=0; u<urgency.size(); u++)
            if(urgency[u] < EPSILON or urgency[u] >= DBL_MAX - EPSILON){
              customerList.erase(customerList.begin() + u);
              urgency.erase(urgency.begin() + u);
              u--;
            }

          this->sort(customerList, urgency);

          /////ADDED
//          customerList.insert(customerList.begin(), 1);
          /////

          succ = customerList.front();


          double breakFlag = false;
          vector<Operation> operations;
          while(customerList.size() > 0){

              if(refuelFlags[shift.getTrailer()] == true)
                  succ = 1;

              allowedTrailerFlag = false;
              for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                  allowedTrailerFlag = true;

              while(not(
                //max driving duration constraint
                (cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]
                        <= this->drivers[lastTimeWindows.front().second].getMaxDrivingDuration()
                /*or cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0)
                        <= this->drivers[lastTimeWindows.front().second].getMaxDrivingDuration()*/)

                and
                //shift within time window
                 (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]
                        <= lastTimeWindows.front().first.second
                 /*or lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0)
                        <= lastTimeWindows.front().first.second*/)

                and
                //shifts with same driver should be distanciated by a minimum intershift duration
                ( (lastTimeWindows[twD].first.first - (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0])
                     >= drivers[lastTimeWindows.front().second].getMinInterSHIFTDURATION()
                  or lastDShift   )
                  /*or (lastTimeWindows[twD].first.first - (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0))
                      >= drivers[lastTimeWindows.front().second].getMinInterSHIFTDURATION()
                   or lastDShift) */)


                //shift using the same trailer cannot overlap
                and (lastTimeWindows[twT].first.second <= lastTimeWindows.front().first.first or lastTimeWindows[twT].first.first >= lastTimeWindows.front().first.second
                 or lastTShift    )
                //customer site should be accessible by the shift's trailer
                and (allowedTrailerFlag or succ==1)


                )


                or (tankQuantities[succ] + (this->actualQuantity[succ][time + cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity())
                                            >= this->customers[succ].getCapacity()-EPSILON and succ != 1)
                //demand for the customer is completely satisfied
                or (this->horizons[succ].back() + deliveredQuantities[succ] >= -EPSILON and succ != 1)

                     or succ == noCustomer
             ){


                  if(customerList.size() <= 1 or succ == 1){
                      breakFlag = true;
                      break;
                  }

                  customerList.erase(customerList.begin());
                  succ = customerList.front();


                  allowedTrailerFlag = false;
                  for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                    if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                      allowedTrailerFlag = true;
              }

              if(breakFlag)
                  break;


              Operation operation;
              operation.setPoint(succ);
              operation.setArrival(time + cumTime + this->timeMatrices[prec][succ]);

              double quantity, actualTankQuantity = 0.0;

              if(succ == 1){
                  quantity = trailerQuantities[shift.getTrailer()] - this->trailers[shift.getTrailer()].getCapacity();
                  refuelFlags[shift.getTrailer()] = false;
              }
              else{
                  actualTankQuantity = tankQuantities[succ] + (this->actualQuantity[succ][operation.getArrival()] - this->customers[succ].getInitialTankQuantity());

                  if(trailerQuantities[shift.getTrailer()] <= this->customers[succ].getCapacity() - actualTankQuantity)
                      quantity = trailerQuantities[shift.getTrailer()];
                  else
                      quantity = this->customers[succ].getCapacity() - actualTankQuantity;

                  if(quantity > abs(this->horizons[succ].back() + deliveredQuantities[succ]))
                      quantity = abs(this->horizons[succ].back() + deliveredQuantities[succ]);
              }

              if(quantity > this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * servingRatio)
                  quantity = this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * (servingRatio);

              operation.setQuantity(quantity);
              trailerQuantities[shift.getTrailer()] -= quantity;
              tankQuantities[succ] += quantity;
              deliveredQuantities[succ] += quantity;

              cumTime += this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime();

              prec = succ;

              customerList.erase(customerList.begin());


//              COUT<<"   "<<operation.getPoint()<<" "<<operation.getArrival()<<" "<<operation.getQuantity()<<"\n";

              vector<unsigned int> customerList;
              for(int c=0; c<this->customers.size(); c++)
                  customerList.push_back(c);

              vector<double> urgency(this->customers.size(), DBL_MAX);

              this->computeUrgency2(customerList, urgency,
                                             shift, deliveredQuantities,
                                             costWeight, demandWeight, prec);



              for(int u=0; u<urgency.size(); u++){
                if(urgency[u] < EPSILON or urgency[u] >= DBL_MAX - EPSILON){
                  customerList.erase(customerList.begin() + u);
                  urgency.erase(urgency.begin() + u);
                  u--;
                }
                if(customerList[u]==succ){
                    customerList.erase(customerList.begin() + u);
                    urgency.erase(urgency.begin() + u);
                    u--;
                }
              }

              COUT<<"\n";
              for(int i=0; i<urgency.size(); i++)
                  COUT<<urgency[i]<<" ";
              COUT<<"\n";
              this->sort(customerList, urgency);

              succ = customerList.front();

              if(trailerQuantities[shift.getTrailer()] < EPSILON + this->trailers[shift.getTrailer()].getCapacity() * refuelRatio)
                  refuelFlags[shift.getTrailer()] = true;

              operations.push_back(operation);
          }

          if(operations.size() != 0){
            shift.setOperations(operations);
            shifts.push_back(shift);
            index++;
   //         COUT<<"SI: "<<shift.getIndex()<<"\n";
            if(shift.getIndex()+1 >= maxShifts)
                break;
          }
          lastTimeWindows.erase(lastTimeWindows.begin());
    }

    solution.setShifts(shifts);


    vector<Shift> fullShifts = solution.getShifts();
    for(int s=fullShifts.size() - 1; s>= 0; s--){
        if(fullShifts[s].getOperations().size() == 0){
            fullShifts.erase(fullShifts.end());
        }

    }
    solution.setShifts(fullShifts);
    fullShifts = solution.getShifts();
/*
        for(int s=0; s<fullShifts.size(); s++){
            COUT<<"SHIFT C: "<<fullShifts[s].getIndex()<<" "<<fullShifts[s].getStart()<<" "<<fullShifts[s].getDriver()<<" "<<fullShifts[s].getTrailer()<<"\n";
            for(int o=0; o<fullShifts[s].getOperations().size(); o++)
              COUT<<"   "<<fullShifts[s].getOperations()[o].getArrival()<<" "
                 <<fullShifts[s].getOperations()[o].getPoint()<<" "
                <<fullShifts[s].getOperations()[o].getQuantity()<<"\n";
        }
*/
    return solution;
}

void Instance::computeUrgency(vector<unsigned int> &customerList, vector<double> &urgency,
                               Shift shift, vector<double> deliveredQuantities,
                               double timeWeight, double quantityWeight, int ties,
                               unsigned int prec){
    for(int c=0; c<this->customers.size(); c++){
        if(c >= 1){
          Customer customer = this->customers[c];
          vector<double> forecast = customer.getForecast();
          double totForecast = 0;
          bool flag = true;

          double deliveredQuantity = deliveredQuantities[c];
          unsigned int currentPoint = prec;

          double totalDistance = 0.0;
          totalDistance = *max_element(this->timeMatrices[currentPoint].begin(), this->timeMatrices[currentPoint].end());
          double totalForecast = this->actualQuantity[c].back() - this->customers[c].getInitialTankQuantity();

          unsigned int f = 0;
          double timeFactor, quantityFactor, distanceFactor;
          while(f < forecast.size() and flag){
            totForecast += forecast[f];
            if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
                timeFactor = ((double)f/forecast.size()) * timeWeight/2;
                quantityFactor = (abs(totalForecast - deliveredQuantity)/totalForecast) * quantityWeight/2;
                distanceFactor = 100.0 * EPSILON * ((double)this->timeMatrices[currentPoint][c]/(totalDistance))*ties;
                urgency[c] = 1.0 + timeFactor + quantityFactor + distanceFactor;
                flag = false;
          }
            f++;
          }
        }
      else
          urgency[c] == DBL_MAX;
    }
}

irpSolution Instance::constructSolution(irpSolution initialSolution, double timeWeight, double quantityWeight, int ties,
                                        double servingRatio, double refuelRatio,
                                        unsigned int noCustomer,
                                        unsigned int maxShifts){

    irpSolution solution = initialSolution;

    vector<Shift> shifts = solution.getShifts();

    if(shifts.size() > maxShifts)
        return solution;

    vector<double> trailerQuantities;
    for(int t=0; t<this->trailers.size(); t++)
        trailerQuantities.push_back(this->trailers[t].getInitialQuantity());

    vector<double> tankQuantities;
    for(int c=0; c<this->customers.size(); c++)
        tankQuantities.push_back(this->customers[c].getInitialTankQuantity());


    vector<double> deliveredQuantities(this->customers.size(), 0.0);

    vector<Operation> lastOperations(this->customers.size(), Operation());
    vector<unsigned int> lastOperationsTime(this->customers.size(), 0);
    unsigned int lastTime = 0;
    for(int s=0; s<shifts.size(); s++){
        Shift shift = shifts[s];
        for(int o=0; o<shift.getOperations().size(); o++){
            Operation operation = shift.getOperations()[o];
            trailerQuantities[shift.getTrailer()] -= operation.getQuantity();
            tankQuantities[operation.getPoint()] += operation.getQuantity();
            deliveredQuantities[operation.getPoint()] += operation.getQuantity();

            if(operation.getArrival() > lastOperationsTime[operation.getPoint()]){
                lastOperationsTime[operation.getPoint()] = operation.getArrival();
                lastOperations[operation.getPoint()] = operation;
            }

        }
        if(shift.getStart()  > lastTime)
            lastTime = shift.getStart();
    }

    vector< pair< pair<unsigned int,unsigned int>, unsigned int > > lastTimeWindows;

    for(int t=0; t<this->timeWindows.size(); t++){
      if(this->timeWindows[t].first.first >= lastTime){

          pair< pair<unsigned int,unsigned int>, unsigned int > tw;
          tw.first = this->timeWindows[t].first;
          tw.second = this->timeWindows[t].second;
          lastTimeWindows.push_back(tw);

          int numberOfDrivers = this->drivers.size();

          if(numberOfDrivers > shifts.size())
              numberOfDrivers -= (numberOfDrivers - shifts.size());

          for(int d=0; d<numberOfDrivers; d++){
              Operation lastOperation = shifts[/*s*/shifts.size()-1-d].getOperations().back();
              int lastOperationTime = lastOperation.getArrival()
                      + this->customers[lastOperation.getPoint()].getSetupTime()
                      + this->timeMatrices[lastOperation.getPoint()][0];
              if( (abs(lastOperationTime - (int)tw.first.first) < this->drivers[tw.second].getMinInterSHIFTDURATION() or shifts[/*s*/initialSolution.getShifts().size()-1-d].getStart()==tw.first.first)
                     and (int)tw.second == (int)shifts[/*s*/shifts.size()-1-d].getDriver() )
                  lastTimeWindows.erase(lastTimeWindows.end());
          }
      }
    }

    if(lastTimeWindows.size() == 0)
        return solution;

    unsigned int time, index = 0;
    if(shifts.size() != 0){
        time = lastTimeWindows.front().first.first;
        index = shifts.back().getIndex() + 1;
    }


    vector<bool> refuelFlags(this->trailers.size(), false);
    for(int t=0; t<this->trailers.size(); t++)
        if(trailerQuantities[t] < EPSILON)
            refuelFlags[t] = true;
//    vector<Shift> shifts = solution.getShifts();

    while(lastTimeWindows.size() > 0)  {

        Shift shift;
        shift.setIndex(index);
        shift.setDriver(lastTimeWindows.front().second);
        shift.setTrailer(this->drivers[shift.getDriver()].getTrailer());
        shift.setStart(lastTimeWindows.front().first.first);

//        COUT<<"SHIFT C: "<<shift.getIndex()<<" "<<shift.getStart()<<" "<<shift.getDriver()<<" "<<shift.getTrailer()<<"\n";

        time = shift.getStart();

        unsigned int prec;
        unsigned int succ;

        prec = 0;
        succ = 1;

        unsigned int cumTime = 0;

        unsigned int twD=0;
        bool lastDShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(lastTimeWindows[t].second == shift.getDriver()){
            twD = t;
            lastDShift = false;
            break;
          }
        unsigned int twT=0;
        bool lastTShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(this->drivers[lastTimeWindows[t].second].getTrailer() == this->drivers[shift.getDriver()].getTrailer()){
            twT = t;
            lastTShift = false;
            break;
          }
          bool allowedTrailerFlag = false;
          for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
            if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[shift.getDriver()].getTrailer())
                allowedTrailerFlag = true;


          vector<unsigned int> customerList;
          for(int c=0; c<this->customers.size(); c++)
              customerList.push_back(c);
          //////
          vector<double> urgency(this->customers.size(),0.0);


          //////////////////////////
          this->computeUrgency(customerList, urgency,
                                         shift, deliveredQuantities,
                                         timeWeight, quantityWeight, ties,
                                         0.0);
          //////////////////////////


          for(int u=0; u<urgency.size(); u++)
            if(urgency[u] <= EPSILON or urgency[u] > DBL_MAX - EPSILON){
              customerList.erase(customerList.begin() + u);
              urgency.erase(urgency.begin() + u);
              u--;
            }

          this->sort(customerList, urgency);
          /////ADDED
//          customerList.insert(customerList.begin(), 1);
          /////

          succ = customerList.front();


          double breakFlag = false;
          vector<Operation> operations;
          while(customerList.size() > 0){

              if(refuelFlags[shift.getTrailer()] == true)
                  succ = 1;

              allowedTrailerFlag = false;
              for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                  allowedTrailerFlag = true;

              while(not(
                //max driving duration constraint
                (cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]
                        <= this->drivers[lastTimeWindows.front().second].getMaxDrivingDuration()
                /*or cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0)
                        <= this->drivers[lastTimeWindows.front().second].getMaxDrivingDuration()*/)

                and
                //shift within time window
                 (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]
                        <= lastTimeWindows.front().first.second
                 /*or lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0)
                        <= lastTimeWindows.front().first.second*/)

                and
                //shifts with same driver should be distanciated by a minimum intershift duration
                ( (lastTimeWindows[twD].first.first - (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0])
                     >= drivers[lastTimeWindows.front().second].getMinInterSHIFTDURATION()
                  or lastDShift   )
                  /*or (lastTimeWindows[twD].first.first - (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0))
                      >= drivers[lastTimeWindows.front().second].getMinInterSHIFTDURATION()
                   or lastDShift) */)


                //shift using the same trailer cannot overlap
                and (lastTimeWindows[twT].first.second <= lastTimeWindows.front().first.first or lastTimeWindows[twT].first.first >= lastTimeWindows.front().first.second
                 or lastTShift    )
                //customer site should be accessible by the shift's trailer
                and (allowedTrailerFlag or succ==1)


                )


                or (tankQuantities[succ] + (this->actualQuantity[succ][time + cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity())
                                            >= this->customers[succ].getCapacity()-EPSILON and succ != 1)
                //demand for the customer is completely satisfied
                or (this->horizons[succ].back() + deliveredQuantities[succ] >= -EPSILON and succ != 1)

                     or succ == noCustomer
             ){


                  if(customerList.size() <= 1 or succ == 1){
                      breakFlag = true;
                      break;
                  }

                  customerList.erase(customerList.begin());
                  succ = customerList.front();


                  allowedTrailerFlag = false;
                  for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                    if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                      allowedTrailerFlag = true;
              }

              if(breakFlag)
                  break;


              Operation operation;
              operation.setPoint(succ);
              operation.setArrival(time + cumTime + this->timeMatrices[prec][succ]);

              double quantity, actualTankQuantity = 0.0;

              if(succ == 1){
                  quantity = trailerQuantities[shift.getTrailer()] - this->trailers[shift.getTrailer()].getCapacity();
                  refuelFlags[shift.getTrailer()] = false;
              }
              else{
                  actualTankQuantity = tankQuantities[succ] + (this->actualQuantity[succ][operation.getArrival()] - this->customers[succ].getInitialTankQuantity());

                  if(trailerQuantities[shift.getTrailer()] <= this->customers[succ].getCapacity() - actualTankQuantity)
                      quantity = trailerQuantities[shift.getTrailer()];
                  else
                      quantity = this->customers[succ].getCapacity() - actualTankQuantity;

                  if(quantity > abs(this->horizons[succ].back() + deliveredQuantities[succ]))
                      quantity = abs(this->horizons[succ].back() + deliveredQuantities[succ]);
              }

              if(quantity > this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * servingRatio)
                  quantity = this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * (servingRatio);

              operation.setQuantity(quantity);
              trailerQuantities[shift.getTrailer()] -= quantity;
              tankQuantities[succ] += quantity;
              deliveredQuantities[succ] += quantity;

              cumTime += this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime();

              prec = succ;

              customerList.erase(customerList.begin());


//              COUT<<"   "<<operation.getPoint()<<" "<<operation.getArrival()<<" "<<operation.getQuantity()<<"\n";

              vector<unsigned int> customerList;
              for(int c=0; c<this->customers.size(); c++)
                  customerList.push_back(c);

              vector<double> urgency(this->customers.size(),0.0);

              //////////////////////////
              this->computeUrgency(customerList, urgency,
                                             shift, deliveredQuantities,
                                             timeWeight, quantityWeight, ties,
                                             prec);
              //////////////////////////
              for(int u=0; u<urgency.size(); u++){
                if(urgency[u] <= EPSILON or urgency[u] > DBL_MAX - EPSILON){
                  customerList.erase(customerList.begin() + u);
                  urgency.erase(urgency.begin() + u);
                  u--;
                }
                if(customerList[u]==succ){
                    customerList.erase(customerList.begin() + u);
                    urgency.erase(urgency.begin() + u);
                    u--;
                }
              }

              this->sort(customerList, urgency);

              succ = customerList.front();

              if(trailerQuantities[shift.getTrailer()] < EPSILON + this->trailers[shift.getTrailer()].getCapacity() * refuelRatio)
                  refuelFlags[shift.getTrailer()] = true;

              operations.push_back(operation);
          }

          if(operations.size() != 0){
            shift.setOperations(operations);
            shifts.push_back(shift);
            index++;
   //         COUT<<"SI: "<<shift.getIndex()<<"\n";
            if(shift.getIndex()+1 >= maxShifts)
                break;
          }
          lastTimeWindows.erase(lastTimeWindows.begin());
    }

    solution.setShifts(shifts);


    vector<Shift> fullShifts = solution.getShifts();
    for(int s=fullShifts.size() - 1; s>= 0; s--){
        if(fullShifts[s].getOperations().size() == 0){
            fullShifts.erase(fullShifts.end());
        }

    }
    solution.setShifts(fullShifts);
    fullShifts = solution.getShifts();
/*
        for(int s=0; s<fullShifts.size(); s++){
            COUT<<"SHIFT C: "<<fullShifts[s].getIndex()<<" "<<fullShifts[s].getStart()<<" "<<fullShifts[s].getDriver()<<" "<<fullShifts[s].getTrailer()<<"\n";
            for(int o=0; o<fullShifts[s].getOperations().size(); o++)
              COUT<<"   "<<fullShifts[s].getOperations()[o].getArrival()<<" "
                 <<fullShifts[s].getOperations()[o].getPoint()<<" "
                <<fullShifts[s].getOperations()[o].getQuantity()<<"\n";
        }
*/
    return solution;
}



irpSolution Instance::randomizedConstructSolution(irpSolution initialSolution, double timeWeight, double quantityWeight, int ties,
                                        double servingRatio, double refuelRatio,
                                        unsigned int noCustomer,
                                        unsigned int maxShifts,
                                        double randomPick,
                                        unsigned int urgencyPolicy){

    irpSolution solution = initialSolution;

    vector<Shift> shifts = solution.getShifts();

    if(shifts.size() > maxShifts)
        return solution;

    vector<double> trailerQuantities;
    for(int t=0; t<this->trailers.size(); t++)
        trailerQuantities.push_back(this->trailers[t].getInitialQuantity());

    vector<double> tankQuantities;
    for(int c=0; c<this->customers.size(); c++)
        tankQuantities.push_back(this->customers[c].getInitialTankQuantity());


    vector<double> deliveredQuantities(this->customers.size(), 0.0);

    vector<Operation> lastOperations(this->customers.size(), Operation());
    vector<unsigned int> lastOperationsTime(this->customers.size(), 0);
    unsigned int lastTime = 0;
    for(int s=0; s<shifts.size(); s++){
        Shift shift = shifts[s];
        for(int o=0; o<shift.getOperations().size(); o++){
            Operation operation = shift.getOperations()[o];
            trailerQuantities[shift.getTrailer()] -= operation.getQuantity();
            tankQuantities[operation.getPoint()] += operation.getQuantity();
            deliveredQuantities[operation.getPoint()] += operation.getQuantity();

            if(operation.getArrival() > lastOperationsTime[operation.getPoint()]){
                lastOperationsTime[operation.getPoint()] = operation.getArrival();
                lastOperations[operation.getPoint()] = operation;
            }

        }
        if(shift.getStart()  > lastTime)
            lastTime = shift.getStart();
    }

    vector< pair< pair<unsigned int,unsigned int>, unsigned int > > lastTimeWindows;

    for(int t=0; t<this->timeWindows.size(); t++){
      if(this->timeWindows[t].first.first >= lastTime){

          pair< pair<unsigned int,unsigned int>, unsigned int > tw;
          tw.first = this->timeWindows[t].first;
          tw.second = this->timeWindows[t].second;
          lastTimeWindows.push_back(tw);

          int numberOfDrivers = this->drivers.size();

          if(numberOfDrivers > shifts.size())
              numberOfDrivers -= (numberOfDrivers - shifts.size());

          for(int d=0; d<numberOfDrivers; d++){
              Operation lastOperation = shifts[/*s*/shifts.size()-1-d].getOperations().back();
              int lastOperationTime = lastOperation.getArrival()
                      + this->customers[lastOperation.getPoint()].getSetupTime()
                      + this->timeMatrices[lastOperation.getPoint()][0];
              if( (abs(lastOperationTime - (int)tw.first.first) < this->drivers[tw.second].getMinInterSHIFTDURATION() or initialSolution.getShifts()[/*s*/shifts.size()-1-d].getStart()==tw.first.first)
                     and (int)tw.second == (int)shifts[/*s*/shifts.size()-1-d].getDriver() )
                  lastTimeWindows.erase(lastTimeWindows.end());
          }
      }
    }

    if(lastTimeWindows.size() == 0)
        return solution;

    unsigned int time, index = 0;
    if(shifts.size() != 0){
        time = lastTimeWindows.front().first.first;
        index = shifts.back().getIndex() + 1;
    }


    vector<bool> refuelFlags(this->trailers.size(), false);
    for(int t=0; t<this->trailers.size(); t++)
        if(trailerQuantities[t] < EPSILON)
            refuelFlags[t] = true;
//    vector<Shift> shifts = shifts();

    while(lastTimeWindows.size() > 0)  {

        Shift shift;
        shift.setIndex(index);
        shift.setDriver(lastTimeWindows.front().second);
        shift.setTrailer(this->drivers[shift.getDriver()].getTrailer());
        shift.setStart(lastTimeWindows.front().first.first);

//        COUT<<"SHIFT C: "<<shift.getIndex()<<" "<<shift.getStart()<<" "<<shift.getDriver()<<" "<<shift.getTrailer()<<"\n";

        time = shift.getStart();

        unsigned int prec;
        unsigned int succ;

        prec = 0;
        succ = 1;

        unsigned int cumTime = 0;

        unsigned int twD=0;
        bool lastDShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(lastTimeWindows[t].second == shift.getDriver()){
            twD = t;
            lastDShift = false;
            break;
          }
        unsigned int twT=0;
        bool lastTShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(this->drivers[lastTimeWindows[t].second].getTrailer() == this->drivers[shift.getDriver()].getTrailer()){
            twT = t;
            lastTShift = false;
            break;
          }
          bool allowedTrailerFlag = false;
          for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
            if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[shift.getDriver()].getTrailer())
                allowedTrailerFlag = true;


          vector<unsigned int> customerList;
          for(int c=0; c<this->customers.size(); c++)
              customerList.push_back(c);
          //////
          vector<double> urgency(this->customers.size(),0.0);
/*
          for(int c=0; c<this->customers.size(); c++){
              if(c >= 1){
                Customer customer = this->customers[c];
                vector<double> forecast = customer.getForecast();
                double totForecast = 0;
                bool flag = true;

                double deliveredQuantity = deliveredQuantities[c];
                unsigned int currentPoint = 0;

                double totalDistance = 0.0;
                totalDistance = *max_element(this->timeMatrices[currentPoint].begin(), this->timeMatrices[currentPoint].end());
                double totalForecast = this->actualQuantity[c].back() - this->customers[c].getInitialTankQuantity();

                unsigned int f = 0;
                double timeFactor, quantityFactor, distanceFactor;
                while(f < forecast.size() and flag){
                  totForecast += forecast[f];
                  if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
                      timeFactor = ((double)f/forecast.size()) * timeWeight/2;
                      quantityFactor = (abs(totalForecast - deliveredQuantity)/totalForecast) * quantityWeight/2;
                      distanceFactor = 100.0 * EPSILON * ((double)this->timeMatrices[currentPoint][c]/(totalDistance))*ties;
                      urgency[c] = 1.0 + timeFactor + quantityFactor + distanceFactor;
                      flag = false;
                }
                  f++;
                }
              }
            else
                urgency[c] == DBL_MAX;
          }
*/
          if(urgencyPolicy == 1)
            this->computeUrgency(customerList, urgency,
                                         shift, deliveredQuantities,
                                         timeWeight, quantityWeight, ties,
                                         0);
          else
            this->computeUrgency2(customerList, urgency,
                                         shift, deliveredQuantities,
                                         timeWeight, quantityWeight,
                                  0);

          for(int u=0; u<urgency.size(); u++)
            if(urgency[u] < EPSILON or urgency[u] >= DBL_MAX - EPSILON){
              customerList.erase(customerList.begin() + u);
              urgency.erase(urgency.begin() + u);
              u--;
            }

          this->sort(customerList, urgency);
          /////ADDED
//          customerList.insert(customerList.begin(), 1);
          /////

          succ = customerList.front();


          double breakFlag = false;
          vector<Operation> operations;
          while(customerList.size() > 0){

              if(refuelFlags[shift.getTrailer()] == true)
                  succ = 1;

              allowedTrailerFlag = false;
              for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                  allowedTrailerFlag = true;

              double rndPick, randomNumber;
              rndPick = randomPick;
              randomNumber = emili::generateRealRandomNumber()/*(double)rand()/RAND_MAX*/;
   //           COUT<<"RANDOM: "<<rndPick<<" "<<randomNumber<<"\n"; //int a; CIN>>a;

              while(not(
                //max driving duration constraint
                (cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]
                        <= this->drivers[lastTimeWindows.front().second].getMaxDrivingDuration()
                /*or cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0)
                        <= this->drivers[lastTimeWindows.front().second].getMaxDrivingDuration()*/)

                and
                //shift within time window
                 (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]
                        <= lastTimeWindows.front().first.second
                 /*or lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0)
                        <= lastTimeWindows.front().first.second*/)

                and
                //shifts with same driver should be distanciated by a minimum intershift duration
                ( (lastTimeWindows[twD].first.first - (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0])
                     >= drivers[lastTimeWindows.front().second].getMinInterSHIFTDURATION()
                  or lastDShift   )
                  /*or (lastTimeWindows[twD].first.first - (lastTimeWindows.front().first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->computeTime(succ, 0))
                      >= drivers[lastTimeWindows.front().second].getMinInterSHIFTDURATION()
                   or lastDShift) */)


                //shift using the same trailer cannot overlap
                and (lastTimeWindows[twT].first.second <= lastTimeWindows.front().first.first or lastTimeWindows[twT].first.first >= lastTimeWindows.front().first.second
                 or lastTShift    )
                //customer site should be accessible by the shift's trailer
                and (allowedTrailerFlag or succ==1)


                )


                or (tankQuantities[succ] + (this->actualQuantity[succ][time + cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity())
                                            >= this->customers[succ].getCapacity()-EPSILON and succ != 1)
                //demand for the customer is completely satisfied
                or (this->horizons[succ].back() + deliveredQuantities[succ] >= -EPSILON and succ != 1)

                     or succ == noCustomer
                     or (rndPick < randomNumber and succ != 1)
             ){

                  randomNumber = emili::generateRealRandomNumber()/*(double)rand()/RAND_MAX*/;
                  rndPick += rndPick/2;
//                  COUT<<"RANDOM: "<<rndPick<<" "<<randomNumber<<"\n";// int a; CIN>>a;

                  if(customerList.size() <= 1 or succ == 1){
                      breakFlag = true;
                      break;
                  }

                  customerList.erase(customerList.begin());
                  succ = customerList.front();


                  allowedTrailerFlag = false;
                  for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                    if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                      allowedTrailerFlag = true;
              }

              if(breakFlag)
                  break;


              Operation operation;
              operation.setPoint(succ);
              operation.setArrival(time + cumTime + this->timeMatrices[prec][succ]);

              double quantity, actualTankQuantity = 0.0;

              if(succ == 1){
                  quantity = trailerQuantities[shift.getTrailer()] - this->trailers[shift.getTrailer()].getCapacity();
                  refuelFlags[shift.getTrailer()] = false;
              }
              else{
                  actualTankQuantity = tankQuantities[succ] + (this->actualQuantity[succ][operation.getArrival()] - this->customers[succ].getInitialTankQuantity());

                  if(trailerQuantities[shift.getTrailer()] <= this->customers[succ].getCapacity() - actualTankQuantity)
                      quantity = trailerQuantities[shift.getTrailer()];
                  else
                      quantity = this->customers[succ].getCapacity() - actualTankQuantity;

                  if(quantity > abs(this->horizons[succ].back() + deliveredQuantities[succ]))
                      quantity = abs(this->horizons[succ].back() + deliveredQuantities[succ]);
              }

              if(quantity > this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * servingRatio)
                  quantity = this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * (servingRatio);

              operation.setQuantity(quantity);
              trailerQuantities[shift.getTrailer()] -= quantity;
              tankQuantities[succ] += quantity;
              deliveredQuantities[succ] += quantity;

              cumTime += this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime();

              prec = succ;

              customerList.erase(customerList.begin());


//              COUT<<"   "<<operation.getPoint()<<" "<<operation.getArrival()<<" "<<operation.getQuantity()<<"\n";

              vector<unsigned int> customerList;
              for(int c=0; c<this->customers.size(); c++)
                  customerList.push_back(c);

              vector<double> urgency(this->customers.size(),0.0);
/*
              for(int c=0; c<this->customers.size(); c++){

                  if(c >= 1){
                    Customer customer = this->customers[c];
                    vector<double> forecast = customer.getForecast();
                    double totForecast = 0;
                    bool flag = true;

                    double deliveredQuantity = deliveredQuantities[c];
                    unsigned int currentPoint = prec;

                    double totalDistance = 0.0;
                    totalDistance = *max_element(this->timeMatrices[currentPoint].begin(), this->timeMatrices[currentPoint].end());
                    double totalForecast = this->actualQuantity[c].back() - this->customers[c].getInitialTankQuantity();

                    unsigned int f = 0;
                    double timeFactor, quantityFactor, distanceFactor = 0.0;
                    while(f < forecast.size() and flag){
                      totForecast += forecast[f];
                      if(totForecast > deliveredQuantity + this->customers[c].getInitialTankQuantity() - this->customers[c].getSafetyLevel()){
                          timeFactor = ((double)f/forecast.size()) * timeWeight/2;
                          quantityFactor = (abs(totalForecast - deliveredQuantity)/totalForecast) * quantityWeight/2;
                          distanceFactor = 100.0 * EPSILON * ((double)this->timeMatrices[currentPoint][c]/(totalDistance))*ties;
                          urgency[c] = 1.0 + timeFactor + quantityFactor + distanceFactor;
                          flag = false;
                      }
                      f++;
                    }
                  }
                else
                    urgency[c] == DBL_MAX;
              }
*/
              if(urgencyPolicy == 1)
                this->computeUrgency(customerList, urgency,
                                             shift, deliveredQuantities,
                                             timeWeight, quantityWeight, ties,
                                             prec);
              else
                this->computeUrgency2(customerList, urgency,
                                             shift, deliveredQuantities,
                                             timeWeight, quantityWeight, prec);

              for(int u=0; u<urgency.size(); u++){
                if(urgency[u] < EPSILON or urgency[u] >= DBL_MAX - EPSILON){
                  customerList.erase(customerList.begin() + u);
                  urgency.erase(urgency.begin() + u);
                  u--;
                }
                if(customerList[u]==succ){
                    customerList.erase(customerList.begin() + u);
                    urgency.erase(urgency.begin() + u);
                    u--;
                }
              }

              this->sort(customerList, urgency);

              succ = customerList.front();

              if(trailerQuantities[shift.getTrailer()] < EPSILON + this->trailers[shift.getTrailer()].getCapacity() * refuelRatio)
                  refuelFlags[shift.getTrailer()] = true;

              operations.push_back(operation);
          }

          if(operations.size() != 0){
            shift.setOperations(operations);
            shifts.push_back(shift);
            index++;
            if(index >= maxShifts)
                break;
          }
          lastTimeWindows.erase(lastTimeWindows.begin());
    }

    solution.setShifts(shifts);


    vector<Shift> fullShifts = solution.getShifts();
    for(int s=fullShifts.size() - 1; s>= 0; s--){
        if(fullShifts[s].getOperations().size() == 0){
            fullShifts.erase(fullShifts.end());
        }

    }
    solution.setShifts(fullShifts);
    fullShifts = solution.getShifts();
/*
        for(int s=0; s<fullShifts.size(); s++){
            COUT<<"SHIFT C: "<<fullShifts[s].getIndex()<<" "<<fullShifts[s].getStart()<<" "<<fullShifts[s].getDriver()<<" "<<fullShifts[s].getTrailer()<<"\n";
            for(int o=0; o<fullShifts[s].getOperations().size(); o++)
              COUT<<"   "<<fullShifts[s].getOperations()[o].getArrival()<<" "
                 <<fullShifts[s].getOperations()[o].getPoint()<<" "
                <<fullShifts[s].getOperations()[o].getQuantity()<<"\n";
        }
*/
    return solution;
}


irpSolution Instance::extendSolution(irpSolution &solution, double servingRatio, double refuelRatio,
                                     unsigned int noCustomer,
                                     unsigned int maxShift){

    vector<double> finalTankQuantities(this->customers.size(), 0.0);
    vector<double> finalTrailerQuantities;
    for(int t=0; t<this->trailers.size(); t++)
        finalTrailerQuantities.push_back(this->trailers[t].getInitialQuantity());
    unsigned int lastTime = 0;
    Operation lastOperation;

    for(int s=0; s<solution.getShifts().size(); s++){
        Shift shift = solution.getShifts()[s];
        for(int o=0; o<shift.getOperations().size(); o++){
            Operation operation = shift.getOperations()[o];
            finalTankQuantities[operation.getPoint()] += operation.getQuantity();
            finalTrailerQuantities[shift.getTrailer()] -= operation.getQuantity();
            if(operation.getArrival() + this->customers[operation.getPoint()].getSetupTime() + this->timeMatrices[operation.getPoint()][0] > lastTime){
                lastTime = operation.getArrival() + this->customers[operation.getPoint()].getSetupTime() + this->timeMatrices[operation.getPoint()][0];
                lastOperation = operation;
            }
        }

    }

    lastTime += this->customers[lastOperation.getPoint()].getSetupTime() + this->timeMatrices[lastOperation.getPoint()][0];
    for(int c=0; c<finalTankQuantities.size(); c++)
        finalTankQuantities[c] += (actualQuantity[c][lastTime]);



    vector< pair< pair<unsigned int,unsigned int>, unsigned int > > lastTimeWindows;
    for(int d=0; d<this->drivers.size(); d++){
      for(int t=0; t<this->drivers[d].getTimeWindows().size(); t++){
          if(this->drivers[d].getTimeWindows()[t].first >= lastTime){
              pair< pair<unsigned int,unsigned int>, unsigned int > tw;
              tw.first = this->drivers[d].getTimeWindows()[t];
              tw.second = d;
              lastTimeWindows.push_back(tw);
          }
      }
    }

    if(lastTimeWindows.size() == 0)
        return solution;

    this->sortPair(lastTimeWindows);


    unsigned int driver = lastTimeWindows.front().second;
    unsigned int trailer = this->drivers[driver].getTrailer();

    vector< vector<double> > closestCustomers(this->customers.size());
    for(int c1=0; c1<closestCustomers.size(); c1++){
        vector<double> cost(this->customers.size());
        for(int c2=0; c2<cost.size(); c2++){
            double distanceCost = this->trailers[trailer].getDistanceCost() * this->distMatrices[c1][c2];
            double timeCost = this->drivers[driver].getTimeCost() * this->timeMatrices[c1][c2];
            if(abs(finalTankQuantities[c2] - this->customers[c2].getCapacity()) > EPSILON)
                cost[c2] =  (distanceCost + timeCost)/abs(finalTankQuantities[c2] - this->customers[c2].getCapacity());
            else
                cost[c2] = DBL_MAX;
            if(c1 == c2)
                cost[c2] = DBL_MAX;
        }
        closestCustomers[c1] = cost;

    }

    vector<double> refuelFlag(this->trailers.size(), false);
    for(int t=0; t<refuelFlag.size(); t++)
        if(finalTrailerQuantities[t] < EPSILON)
            refuelFlag[t] = true;


    vector<Shift> shifts;
    unsigned int lastIndex = solution.getShifts().size();
    while(lastTimeWindows.size() > 0){

        vector<unsigned int> customerList;
        for(int c=0; c<this->customers.size(); c++)
          customerList.push_back(c);
        vector<double> urgency(this->customers.size(),0.0);

        unsigned int prec = 0;

        vector< vector<unsigned int> > cc = this->contributeMatrixes.find(driver)->second;
        urgency = closestCustomers[prec];

        this->sort(customerList, urgency);
  //      customerList = cc[prec];


        for(int i=0; i<customerList.size(); i++){
            if(customerList[i] == 0 or customerList[i] == 1){
                customerList.erase(customerList.begin() + i);
                i--;
            }
        }

//        customerList.insert(customerList.begin(), 1);

        driver = lastTimeWindows.front().second;
        trailer = this->drivers[driver].getTrailer();

        Shift shift;
        shift.setIndex(lastIndex);
        shift.setStart(lastTimeWindows.front().first.first);
        shift.setDriver(driver);
        shift.setTrailer(trailer);

        if(finalTrailerQuantities[trailer] < EPSILON ){
            refuelFlag[trailer] = true;
            if(customerList.front() != 1)
                customerList.insert(customerList.begin(), 1);
        }

        unsigned int succ =  customerList.front();

 //       COUT<<"SHIFT: "<<shift.getIndex()<<" "<<shift.getStart()<<" "<<shift.getDriver()<<" "<<shift.getTrailer()<<"\n";
        lastIndex++;

        unsigned int twD=0;
        bool lastDShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(lastTimeWindows[t].second == lastTimeWindows[0].second){
            twD = t;
            lastDShift = false;
            break;
          }
        unsigned int twT=0;
        bool lastTShift = true;
        for(int t=1; t<lastTimeWindows.size(); t++)
          if(this->drivers[lastTimeWindows[t].second].getTrailer() == this->drivers[lastTimeWindows[0].second].getTrailer()){
            twT = t;
            lastTShift = false;
            break;
          }
          bool allowedTrailerFlag = false;
          for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
            if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
              allowedTrailerFlag = true;

        unsigned int cumTime = 0;
        unsigned int time = lastTimeWindows.front().first.first;
        bool breakFlag = false;


        vector<Operation> operations;
        while(customerList.size() > 0){

            if(finalTrailerQuantities[trailer] < EPSILON ){
                refuelFlag[trailer] = true;
                if(customerList.front() != 1)
                    customerList.insert(customerList.begin(), 1);
            }
            succ =  customerList.front();

            allowedTrailerFlag = false;

            for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
              if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                allowedTrailerFlag = true;


            while(not(
              //max driving duration constraint
              cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0] <= this->drivers[lastTimeWindows.front().second].getMaxDrivingDuration()
              //shift within time window
              and lastTimeWindows[0].first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0] <= lastTimeWindows[0].first.second
              //shifts with same driver should be distanciated by a minimum intershift duration
              and (lastTimeWindows[twD].first.first - (lastTimeWindows[0].first.first + cumTime + this->timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) >= drivers[lastTimeWindows[0].second].getMinInterSHIFTDURATION()
                or lastDShift   )
              //shift using the same trailer cannot overlap
              and (lastTimeWindows[twT].first.second <= lastTimeWindows[0].first.first or lastTimeWindows[twT].first.first >= lastTimeWindows[0].first.second
               or lastTShift    )
              //customer site should be accessible by the shift's trailer
              and (allowedTrailerFlag or succ==1)
              )

              or (finalTankQuantities[succ] >= this->customers[succ].getCapacity()-EPSILON and succ != 1)

                  or succ == noCustomer

           ){

                  //it is not possible to go to the source to refuel, shift ends
                  if(succ == 1){
                      breakFlag = true;
                      break;
                  }
                  customerList.erase(customerList.begin());

                  //all customers have been served or cannot be served, shift ends
                  if(customerList.size() <= 1 ){
                    breakFlag = true;
                    break;
                  }

           //      prec = succ;
                 succ = customerList.front();

                 allowedTrailerFlag = false;
                 for(int t=0; t<this->customers[succ].getAllowedTrailers().size(); t++)
                   if(this->customers[succ].getAllowedTrailers()[t] == this->drivers[lastTimeWindows[0].second].getTrailer()  or succ <= 1)
                     allowedTrailerFlag = true;

            }

            if(breakFlag)
                break;

            cumTime += this->timeMatrices[prec][succ];
            time = lastTimeWindows.front().first.first + cumTime;

            Operation operation;
            operation.setPoint(succ);
            operation.setArrival(time);

//            COUT<<"                     "<<finalTrailerQuantities[trailer]<<" "<<finalTankQuantities[succ]<<" "<<this->customers[succ].getCapacity()<<"\n";
            double quantity = 0.0;
            if(succ == 1){
                quantity = finalTrailerQuantities[trailer] - this->trailers[trailer].getCapacity();
                refuelFlag[trailer] = false;
            }
            else{
                if(finalTrailerQuantities[trailer] <= (this->customers[succ].getCapacity() - finalTankQuantities[succ]))
                    quantity = finalTrailerQuantities[trailer];
                else
                    quantity = (this->customers[succ].getCapacity() - finalTankQuantities[succ]);
            }

            if(quantity > this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * servingRatio)
                quantity = this->trailers[this->drivers[lastTimeWindows[0].second].getTrailer()].getCapacity() * (servingRatio);

            finalTrailerQuantities[trailer] -= quantity;
            finalTankQuantities[succ] += quantity;
            operation.setQuantity(quantity);

            if((quantity >= EPSILON and succ !=1) or (quantity <= -EPSILON and succ == 1 or operations.size() == 0))
                operations.push_back(operation);

            cumTime +=  this->customers[succ].getSetupTime();

            prec = succ;

            customerList = vector<unsigned int>();
            for(int c=0; c<this->customers.size(); c++)
              customerList.push_back(c);

            vector< vector<unsigned int> > cc = this->contributeMatrixes.find(driver)->second;

            urgency = closestCustomers[prec];


            this->sort(customerList, urgency);
 //           customerList = cc[prec];

            for(int i=0; i<customerList.size(); i++){
                if(customerList[i] == 0 or customerList[i] == 1){
                    customerList.erase(customerList.begin() + i);
                    i--;
                }
            }

            if(finalTrailerQuantities[trailer] < EPSILON + this->trailers[trailer].getCapacity() * refuelRatio){
                refuelFlag[trailer] = true;
                if(customerList.front() != 1)
                    customerList.insert(customerList.begin(), 1);
            }

            succ =  customerList.front();
            }


        shift.setOperations(operations);
        shifts.push_back(shift);

        if(lastIndex-1 >= maxShift)
            break;


        lastTimeWindows.erase(lastTimeWindows.begin());

    }
/*
    COUT<<"NEW SHIFTS: \n";
    for(int s=0; s<shifts.size(); s++){
        COUT<<"SHIFT: "<<shifts[s].getIndex()<<" "<<shifts[s].getStart()<<" "<<shifts[s].getDriver()<<" "<<shifts[s].getTrailer()<<"\n";
        for(int o=0; o<shifts[s].getOperations().size(); o++)
          COUT<<"   "<<shifts[s].getOperations()[o].getArrival()<<" "<<shifts[s].getOperations()[o].getPoint()<<" "<<shifts[s].getOperations()[o].getQuantity()<<"\n";
    }

*/
    vector<Shift> allShifts;
    allShifts = solution.getShifts();
    for(int s=0; s<shifts.size(); s++)
        allShifts.push_back(shifts[s]);
    solution.setShifts(allShifts);

    allShifts = solution.getShifts();
    while(allShifts.back().getOperations().size() == 0
          or (allShifts.back().getOperations().size() == 1 and allShifts.back().getOperations().back().getPoint() == 1))
        allShifts.erase(allShifts.end());
    solution.setShifts(allShifts);

/*
    COUT<<"FULL SHIFTS: \n";
    for(int s=0; s<solution.getShifts().size(); s++){
        COUT<<"SHIFT: "<<solution.getShifts()[s].getIndex()<<" "<<solution.getShifts()[s].getStart()<<" "<<solution.getShifts()[s].getDriver()<<" "<<solution.getShifts()[s].getTrailer()<<"\n";
        for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++)
          COUT<<"   "<<solution.getShifts()[s].getOperations()[o].getArrival()<<" "
             <<solution.getShifts()[s].getOperations()[o].getPoint()<<" "
            <<solution.getShifts()[s].getOperations()[o].getQuantity()<<"\n";
    }
*/

    solution.saveSolution("ExtendedSolution.xml");

    return solution;
}

irpSolution Instance::exchangeShiftSolution(irpSolution solution, unsigned int shiftIndex, unsigned int operationIndex, double servingRatio, double refuelRatio){

    if(shiftIndex >= solution.getShifts().size())
        return solution;
    if(operationIndex >= solution.getShifts()[shiftIndex].getOperations().size() - 1)
        return solution;


    vector<Shift> shifts = solution.getShifts();
    Shift shift = shifts[shiftIndex];
    if(shift.getOperations().size() < 2)
        return solution;

    vector<Operation> operations = shift.getOperations();
    unsigned int departurePoint;
    unsigned int departureTime;
    if(operationIndex == 0){
        departurePoint = 0;
        departureTime = shift.getStart();
    }
    else{
        departurePoint = operations[operationIndex - 1].getPoint();
        departureTime = operations[operationIndex - 1].getArrival() + this->customers[departurePoint].getSetupTime();
    }

    unsigned int point1 = operations[operationIndex].getPoint();
    unsigned int point2 = operations[operationIndex+1].getPoint();

//    if(point1==1 /*or point2==1*/)
//        return solution;

    unsigned int firstTime = departureTime + this->timeMatrices[departurePoint][point2];
    unsigned int secondTime = firstTime + this->customers[point2].getSetupTime() + this->timeMatrices[point2][point1];

    unsigned int arrivalPoint;
    unsigned int arrivalTime;
    if(operationIndex == operations.size()-2){
        arrivalPoint == 0;
        arrivalTime = secondTime + this->customers[point1].getSetupTime() + this->timeMatrices[point1][0];
    }
    else{
        arrivalPoint = operations[operationIndex + 2].getPoint();
        arrivalTime = secondTime + this->customers[point1].getSetupTime() + this->timeMatrices[point1][arrivalPoint];
    }

    double trailerQuantity = this->trailers[shifts[shiftIndex].getTrailer()].getInitialQuantity();
    vector<double> tankQuantities;
    for(int t=0; t<this->customers.size(); t++)
        tankQuantities.push_back(this->customers[t].getInitialTankQuantity());

    unsigned int driver = shifts[shiftIndex].getTrailer();
    unsigned int trailer = shifts[shiftIndex].getTrailer();


    for(int s=0; s<=shiftIndex; s++){
        if(s < shiftIndex)
            for(int o=0; o<shifts[s].getOperations().size(); o++){
                if(shifts[s].getTrailer() == trailer)
                    trailerQuantity -= shifts[s].getOperations()[o].getQuantity();
                tankQuantities[shifts[s].getOperations()[o].getPoint()] += shifts[s].getOperations()[o].getQuantity();
            }
        else
            for(int o=0; o<operationIndex; o++){
                if(shifts[s].getTrailer() == trailer)
                    trailerQuantity -= shifts[s].getOperations()[o].getQuantity();
                tankQuantities[shifts[s].getOperations()[o].getPoint()] += shifts[s].getOperations()[o].getQuantity();
            }
    }


    vector<unsigned int> customerList;
    customerList.push_back(point2);
    customerList.push_back(point1);

    for(int o=operationIndex+2; o<shift.getOperations().size(); o++)
        if(shift.getOperations()[o].getPoint() != 1)
            customerList.push_back(shift.getOperations()[o].getPoint());

    Shift newShift;
    newShift.setIndex(shift.getIndex());
    newShift.setDriver(shift.getDriver());
    newShift.setTrailer(shift.getTrailer());
    newShift.setStart(shift.getStart());

    vector<Operation> newOperations;
    for(int o=0; o<operationIndex; o++)
        newOperations.push_back(operations[o]);

    unsigned int prec, succ;

    unsigned int time;
    if(operationIndex == 0){
        time = shift.getStart();
        prec = 0;
    }
    else{
        prec = newOperations.back().getPoint();
        time = newOperations.back().getArrival() + this->customers[prec].getSetupTime();
    }

    unsigned int cumTime = time;

    unsigned int twd = 0;
    bool twdFlag = false;
    for(int t=shiftIndex+1; t<shifts.size(); t++){
        if(shifts[t].getDriver() == shift.getDriver()){
            twd = t;
            twdFlag = true;
        }
        if(twdFlag)
            break;
    }

    unsigned tw;
    for(int t=0; t<this->timeWindows.size(); t++)
        if(this->timeWindows[t].first.first == shift.getStart() and this->timeWindows[t].second == driver)
            tw = t;


    bool refuelFlag = false;


    while(customerList.size() > 0){

        if(trailerQuantity < EPSILON)
            refuelFlag = true;

        succ = customerList.front();

        bool breakFlag = false;

        unsigned int timeWindowDuration;
        if(tw < this->timeWindows.size())
            timeWindowDuration = this->timeWindows[tw].first.second - this->timeWindows[tw].first.first;
        else
            timeWindowDuration = UINT_MAX;

        double actualTankQuantity = (tankQuantities[succ]) +
                (this->actualQuantity[succ][cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity());

        while((cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) - shift.getStart() >
                this->drivers[driver].getMaxDrivingDuration()
              or
              (cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) - shift.getStart() >
                              timeWindowDuration
              or
              (shifts[twd].getStart() - (cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) <
                this->drivers[driver].getMinInterSHIFTDURATION() and twdFlag)

              or (actualTankQuantity > this->customers[succ].getCapacity() - EPSILON and succ!=1)
              or (trailerQuantity < EPSILON and succ!=1)
           ){

            if(succ == 1){
                breakFlag = true;
                break;
            }
            if(customerList.size() > 1){
                customerList.erase(customerList.begin());
                succ = customerList.front();
            }
            else{
                breakFlag = true;
                break;
            }

            actualTankQuantity = (tankQuantities[succ]) +
                            (this->actualQuantity[succ][cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity());
        }

        if(breakFlag)
            break;

        cumTime += this->timeMatrices[prec][succ];
        
        Operation operation;
        operation.setArrival(cumTime);
        operation.setPoint(succ);


        double quantity = 0.0;
        
        /*double */actualTankQuantity = (tankQuantities[succ]) +
                (this->actualQuantity[succ][operation.getArrival()] - this->customers[succ].getInitialTankQuantity());


        if(succ == 1){
            quantity = -(this->trailers[trailer].getCapacity() - trailerQuantity);
        }
        else{
            if(trailerQuantity < this->customers[succ].getCapacity() - actualTankQuantity){
                    quantity = trailerQuantity;
            }
            else{
                    quantity = this->customers[succ].getCapacity() - actualTankQuantity;
            }
        }

        if(quantity > this->trailers[trailer].getCapacity() * servingRatio)
            quantity = this->trailers[trailer].getCapacity() * servingRatio;

        if(quantity < -EPSILON and succ!=1)
            quantity = 0.0;

        operation.setQuantity(quantity);

        trailerQuantity -= quantity;
        tankQuantities[succ] += quantity;

        cumTime += this->customers[succ].getSetupTime();

        newOperations.push_back(operation);

        prec = succ;



        customerList.erase(customerList.begin());

        if(trailerQuantity < EPSILON + this->trailers[trailer].getCapacity() * refuelRatio)
            customerList.insert(customerList.begin(), 1);

    }

    newShift.setOperations(newOperations);

    solution.setShift(newShift, shiftIndex);

    shifts = solution.getShifts();

    for(int s=0; s<shifts.size(); s++){
        COUT<<"SHIFT E: "<<shifts[s].getIndex()<<" "<<shifts[s].getStart()<<" "<<shifts[s].getDriver()<<"\n";
        for(int o=0; o<shifts[s].getOperations().size(); o++)
            COUT<<"     "<<shifts[s].getOperations()[o].getPoint()<<" "
                  <<shifts[s].getOperations()[o].getArrival()<<" "
                    <<shifts[s].getOperations()[o].getQuantity()<<"\n";
    }

    return solution;

}

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

irpSolution Instance::insertShiftSolution(irpSolution solution, unsigned int shiftIndex, unsigned int operationIndex, unsigned int insertedOperation,
                                          double servingRatio, double refuelRatio){

    if(shiftIndex > solution.getShifts().size())
        return solution;
    if(operationIndex > solution.getShifts()[shiftIndex].getOperations().size()+1)
        return solution;


    Shift shift = solution.getShifts()[shiftIndex];
    vector<Operation> operations = shift.getOperations();


    unsigned int point1, point2;
    if(operationIndex == 0){
        point1 = 0;
        point2 = operations[operationIndex].getPoint();
    }
    else if(operationIndex == operations.size()){
        point1 = operations[operationIndex-1].getPoint();
        point2 = 0;
    }
    else{
        point1 = operations[operationIndex-1].getPoint();
        point2 = operations[operationIndex].getPoint();
    }




    vector<unsigned int> routes = this->fastestSubroute[point1][point2];
/*
    for(int i=0; i<routes.size(); i++)
        COUT<<routes[i]<<" ";
    COUT<<"\n";
*/
    vector<long unsigned int> indexes = sort_indexes(routes);
 /*   for(int i=0; i<indexes.size(); i++) {
      COUT << indexes[i] << " "<<i<<endl;
    }*/
    unsigned int i = 0;
    while(indexes[i] == 0 or indexes[i] == 1 or indexes[i] == point1 or indexes[i] == point2){
        i++;
    }
    if(insertedOperation == 0)
        insertedOperation = indexes[i];

    COUT<<"\nFASTEST: "<<point1<<" "<<point2<<" "<<insertedOperation<<"\n";

    bool allowedTrailerFlag = false;
    for(int t=0; t<this->customers[insertedOperation].getAllowedTrailers().size(); t++)
      if(this->customers[insertedOperation].getAllowedTrailers()[t] == this->drivers[solution.getShifts()[shiftIndex].getDriver()].getTrailer())
          allowedTrailerFlag = true;

    if(not allowedTrailerFlag)
        return solution;

    double trailerQuantity = this->trailers[solution.getShifts()[shiftIndex].getTrailer()].getInitialQuantity();
    vector<double> tankQuantities;
    for(int t=0; t<this->customers.size(); t++)
        tankQuantities.push_back(this->customers[t].getInitialTankQuantity());

    unsigned int driver = solution.getShifts()[shiftIndex].getTrailer();
    unsigned int trailer = solution.getShifts()[shiftIndex].getTrailer();


    for(int s=0; s<=shiftIndex; s++){
        if(s < shiftIndex)
            for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++){
                if(solution.getShifts()[s].getTrailer() == trailer)
                    trailerQuantity -= solution.getShifts()[s].getOperations()[o].getQuantity();
                tankQuantities[solution.getShifts()[s].getOperations()[o].getPoint()] += solution.getShifts()[s].getOperations()[o].getQuantity();
            }
        else
            for(int o=0; o<operationIndex; o++){
                if(solution.getShifts()[s].getTrailer() == trailer)
                    trailerQuantity -= solution.getShifts()[s].getOperations()[o].getQuantity();
                tankQuantities[solution.getShifts()[s].getOperations()[o].getPoint()] += solution.getShifts()[s].getOperations()[o].getQuantity();
            }
    }


    vector<unsigned int> customerList;
    customerList.push_back(insertedOperation);

//    COUT<<"CC: ";
    for(int o=operationIndex; o<shift.getOperations().size(); o++)
        if(shift.getOperations()[o].getPoint() != 1){
            customerList.push_back(shift.getOperations()[o].getPoint());
//            COUT<<customerList.back()<<" ";
        }
//    COUT<<"\n";

    Shift newShift;
    newShift.setIndex(shift.getIndex());
    newShift.setDriver(shift.getDriver());
    newShift.setTrailer(shift.getTrailer());
    newShift.setStart(shift.getStart());

    vector<Operation> newOperations;
    for(int o=0; o<operationIndex; o++)
        newOperations.push_back(operations[o]);

    unsigned int prec, succ;

    unsigned int time;
    if(operationIndex == 0){
        time = shift.getStart();
        prec = 0;
    }
    else{
        prec = newOperations.back().getPoint();
        time = newOperations.back().getArrival() + this->customers[prec].getSetupTime();
    }

    unsigned int cumTime = time;

    unsigned int twd = 0;
    bool twdFlag = false;
    for(int t=shiftIndex+1; t<solution.getShifts().size(); t++){
        if(solution.getShifts()[t].getDriver() == shift.getDriver()){
            twd = t;
            twdFlag = true;
        }
        if(twdFlag)
            break;
    }

    unsigned tw;
    for(int t=0; t<this->timeWindows.size(); t++)
        if(this->timeWindows[t].first.first == shift.getStart() and this->timeWindows[t].second == driver)
            tw = t;


    while(customerList.size() > 0){

        if(trailerQuantity < EPSILON + this->trailers[trailer].getCapacity() * refuelRatio and customerList.front() != 1)
            customerList.insert(customerList.begin(), 1);

        succ = customerList.front();

        bool breakFlag = false;

        unsigned int timeWindowDuration;
        if(tw < this->timeWindows.size())
            timeWindowDuration = this->timeWindows[tw].first.second - this->timeWindows[tw].first.first;
        else
            timeWindowDuration = UINT_MAX;

        double actualTankQuantity = (tankQuantities[succ]) +
                (this->actualQuantity[succ][cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity());

        while((cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) - shift.getStart() >
                this->drivers[driver].getMaxDrivingDuration()
              or
              (cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) - shift.getStart() >
                              timeWindowDuration
              or
              (solution.getShifts()[twd].getStart() - (cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) <
                this->drivers[driver].getMinInterSHIFTDURATION() and twdFlag)

              or (actualTankQuantity > this->customers[succ].getCapacity() - EPSILON and succ!=1)
              or (trailerQuantity < EPSILON and succ!=1)

           ){

            if(succ == 1){
                breakFlag = true;
                break;
            }
            if(customerList.size() > 1){
                customerList.erase(customerList.begin());
                succ = customerList.front();
            }
            else{
                breakFlag = true;
                break;
            }

            actualTankQuantity = (tankQuantities[succ]) +
                            (this->actualQuantity[succ][cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity());

        }

        if(breakFlag)
            break;

        cumTime += this->timeMatrices[prec][succ];

        Operation operation;
        operation.setArrival(cumTime);
        operation.setPoint(succ);


        double quantity = 0.0;

/*        double */actualTankQuantity = (tankQuantities[succ]) +
                (this->actualQuantity[succ][operation.getArrival()] - this->customers[succ].getInitialTankQuantity());


        if(succ == 1){
            quantity = -(this->trailers[trailer].getCapacity() - trailerQuantity);
        }
        else{
            if(trailerQuantity < this->customers[succ].getCapacity() - actualTankQuantity){
                    quantity = trailerQuantity;
            }
            else{
                    quantity = this->customers[succ].getCapacity() - actualTankQuantity;
            }
        }


        if(quantity > this->trailers[trailer].getCapacity() * servingRatio)
            quantity = this->trailers[trailer].getCapacity() * servingRatio;

        if(quantity < -EPSILON and succ!=1)
            quantity = 0.0;

        operation.setQuantity(quantity);

        trailerQuantity -= quantity;
        tankQuantities[succ] += quantity;

        cumTime += this->customers[succ].getSetupTime();

        newOperations.push_back(operation);

        prec = succ;



        customerList.erase(customerList.begin());

        if(trailerQuantity < EPSILON + this->trailers[trailer].getCapacity() * refuelRatio)
            customerList.insert(customerList.begin(), 1);

    }

    newShift.setOperations(newOperations);

    solution.setShift(newShift, shiftIndex);

    vector<Shift> shifts = solution.getShifts();
    unsigned int size = shifts.size();
    for(int s=size-1; s>=0; s--){
        Shift newShift = shifts[s];
        vector<Operation> newOperations = shifts[s].getOperations();
        unsigned int osize = shifts[s].getOperations().size();
        for(int o=osize-1; o>=0; o--)
            if(shifts[s].getOperations()[o].getQuantity() < EPSILON && shifts[s].getOperations()[o].getQuantity() > -EPSILON )
                newOperations.erase(newOperations.end());
            else
                break;

        newShift.setOperations(newOperations);
        if(newShift.getOperations().size() > 0)
         shifts[s] = newShift;
        else
            shifts.erase(shifts.end());
    }
    solution.setShifts(shifts);
/*
    for(int s=0; s<solution.getShifts().size(); s++){
        COUT<<"SHIFT E: "<<solution.getShifts()[s].getIndex()<<" "<<solution.getShifts()[s].getStart()<<" "<<solution.getShifts()[s].getDriver()<<"\n";
        for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++)
            COUT<<"     "<<solution.getShifts()[s].getOperations()[o].getPoint()<<" "
                  <<solution.getShifts()[s].getOperations()[o].getArrival()<<" "
                    <<solution.getShifts()[s].getOperations()[o].getQuantity()<<"\n";
    }
*/
    return solution;

}

irpSolution Instance::removeShiftSolution(irpSolution solution, unsigned int shiftIndex, unsigned int operationIndex, double servingRatio, double refuelRatio){

    if(shiftIndex > solution.getShifts().size())
        return solution;
    if(operationIndex >= solution.getShifts()[shiftIndex].getOperations().size())
        return solution;


    Shift shift = solution.getShifts()[shiftIndex];
    vector<Operation> operations = shift.getOperations();


    unsigned int point1, point2;
    if(operationIndex == 0 and operations.size()==1){
        vector<Shift> newShifts = solution.getShifts();
        newShifts.erase(newShifts.begin() + shiftIndex);
        for(int s = shiftIndex; s<newShifts.size(); s++){
            Shift newShift = newShifts[s];
            newShift.setIndex(newShift.getIndex()-1);
            newShifts[s] = newShift;
        }
        solution.setShifts(newShifts);
        return solution;
    }

    if(operationIndex == 0){
        point1 = 0;
        point2 = operations[operationIndex+1].getPoint();
    }
    else if(operationIndex == operations.size()-1){
        point1 = operations[operationIndex-1].getPoint();
        point2 = 0;
    }
    else{
        point1 = operations[operationIndex-1].getPoint();
        point2 = operations[operationIndex+1].getPoint();
    }


    double trailerQuantity = this->trailers[solution.getShifts()[shiftIndex].getTrailer()].getInitialQuantity();
    vector<double> tankQuantities;
    for(int t=0; t<this->customers.size(); t++)
        tankQuantities.push_back(this->customers[t].getInitialTankQuantity());

    unsigned int driver = solution.getShifts()[shiftIndex].getTrailer();
    unsigned int trailer = solution.getShifts()[shiftIndex].getTrailer();


    for(int s=0; s<=shiftIndex; s++){
        if(s < shiftIndex)
            for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++){
                if(solution.getShifts()[s].getTrailer() == trailer)
                    trailerQuantity -= solution.getShifts()[s].getOperations()[o].getQuantity();
                tankQuantities[solution.getShifts()[s].getOperations()[o].getPoint()] += solution.getShifts()[s].getOperations()[o].getQuantity();
            }
        else{
            unsigned int i = 0;
            if((int)operationIndex> 0)
                i = (int)operationIndex;
            for(int o=0; o<i; o++){
                if(solution.getShifts()[s].getTrailer() == trailer)
                    trailerQuantity -= solution.getShifts()[s].getOperations()[o].getQuantity();
                tankQuantities[solution.getShifts()[s].getOperations()[o].getPoint()] += solution.getShifts()[s].getOperations()[o].getQuantity();
            }
        }
    }


    vector<unsigned int> customerList;

    unsigned int i = operations.size()-1;
    if((int)operationIndex+1 < operations.size())
        i = (int)operationIndex+1;
    COUT<<"CC: ";
    for(int o=i; o<shift.getOperations().size(); o++)
        if(shift.getOperations()[o].getPoint() != 1 and o != operationIndex and o >= 0){
            customerList.push_back(shift.getOperations()[o].getPoint());
            COUT<<customerList.back()<<" ";
        }
    COUT<<"\n";

    Shift newShift;
    newShift.setIndex(shift.getIndex());
    newShift.setDriver(shift.getDriver());
    newShift.setTrailer(shift.getTrailer());
    newShift.setStart(shift.getStart());

    COUT<<"\nOperation: \n";
    vector<Operation> newOperations;
    i = 0;
    if((int)operationIndex > 0)
        i = (int)operationIndex;
    for(int o=0; o<i; o++){
        if(o>=0){
            newOperations.push_back(operations[o]);
            COUT<<newOperations.back().getPoint()<<" ";
        }
    }
//    COUT<<"\n";
    unsigned int prec, succ;

    unsigned int time;
    if(operationIndex == 0){
        time = shift.getStart();
        prec = 0;
    }
    else{
        prec = newOperations.back().getPoint();
        time = newOperations.back().getArrival() + this->customers[prec].getSetupTime();
    }

    unsigned int cumTime = time;

    unsigned int twd = 0;
    bool twdFlag = false;
    for(int t=shiftIndex+1; t<solution.getShifts().size(); t++){
        if(solution.getShifts()[t].getDriver() == shift.getDriver()){
            twd = t;
            twdFlag = true;
        }
        if(twdFlag)
            break;
    }

    unsigned tw;
    for(int t=0; t<this->timeWindows.size(); t++)
        if(this->timeWindows[t].first.first == shift.getStart() and this->timeWindows[t].second == driver)
            tw = t;


    while(customerList.size() > 0){

        if(trailerQuantity < EPSILON + this->trailers[trailer].getCapacity() * refuelRatio and customerList.front()!=1)
            customerList.insert(customerList.begin(), 1);

        succ = customerList.front();

        bool breakFlag = false;

        unsigned int timeWindowDuration;
        if(tw < this->timeWindows.size())
            timeWindowDuration = this->timeWindows[tw].first.second - this->timeWindows[tw].first.first;
        else
            timeWindowDuration = UINT_MAX;


        double actualTankQuantity = (tankQuantities[succ]) +
                (this->actualQuantity[succ][cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity());


        while((cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) - shift.getStart() >
                this->drivers[driver].getMaxDrivingDuration()
              or
              (cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) - shift.getStart() >
                              timeWindowDuration
              or
              (solution.getShifts()[twd].getStart() - (cumTime + timeMatrices[prec][succ] + this->customers[succ].getSetupTime() + this->timeMatrices[succ][0]) <
                this->drivers[driver].getMinInterSHIFTDURATION() and twdFlag)

              or (actualTankQuantity > this->customers[succ].getCapacity() - EPSILON and succ!=1)
              or (trailerQuantity < EPSILON and succ!=1)
           ){


            if(succ == 1){
                breakFlag = true;
                break;
            }
            if(customerList.size() > 1){
                customerList.erase(customerList.begin());
                succ = customerList.front();
            }
            else{
                breakFlag = true;
                break;
            }
            actualTankQuantity = (tankQuantities[succ]) +
                           (this->actualQuantity[succ][cumTime + this->timeMatrices[prec][succ]] - this->customers[succ].getInitialTankQuantity());

        }

        if(breakFlag)
            break;

        cumTime += this->timeMatrices[prec][succ];

        Operation operation;
        operation.setArrival(cumTime);
        operation.setPoint(succ);


        double quantity = 0.0;

        /*double*/ actualTankQuantity = (tankQuantities[succ]) +
                (this->actualQuantity[succ][operation.getArrival()] - this->customers[succ].getInitialTankQuantity());

        if(succ == 1){
            quantity = -(this->trailers[trailer].getCapacity() - trailerQuantity);
        }
        else{
            if(trailerQuantity < this->customers[succ].getCapacity() - actualTankQuantity){
                    quantity = trailerQuantity;
            }
            else{
                    quantity = this->customers[succ].getCapacity() - actualTankQuantity;
            }
        }


        if(quantity > this->trailers[trailer].getCapacity() * servingRatio){
            quantity = this->trailers[trailer].getCapacity() * servingRatio;
            COUT<<operation.getPoint()<<" "<<operation.getArrival()<<" "<<quantity<<"\n";
        }

        if(quantity < -EPSILON and succ!=1){
            quantity = 0.0;
        }

        operation.setQuantity(quantity);

        trailerQuantity -= quantity;
        tankQuantities[succ] += quantity;

        cumTime += this->customers[succ].getSetupTime();

        newOperations.push_back(operation);

        prec = succ;


        customerList.erase(customerList.begin());

        if(trailerQuantity < EPSILON + this->trailers[trailer].getCapacity() * refuelRatio)
            customerList.insert(customerList.begin(), 1);

    }

    newShift.setOperations(newOperations);

    if(newShift.getOperations().size() > 0)
        solution.setShift(newShift, shiftIndex);
/*
    for(int s=0; s<solution.getShifts().size(); s++){
        COUT<<"SHIFT E: "<<solution.getShifts()[s].getIndex()<<" "<<solution.getShifts()[s].getStart()<<" "<<solution.getShifts()[s].getDriver()<<"\n";
        for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++)
            COUT<<"     "<<solution.getShifts()[s].getOperations()[o].getPoint()<<" "
                  <<solution.getShifts()[s].getOperations()[o].getArrival()<<" "
                    <<solution.getShifts()[s].getOperations()[o].getQuantity()<<"\n";
    }
*/
    return solution;

}

irpSolution Instance::perturbation(irpSolution solution, unsigned int customer){

    irpSolution newSolution = solution;

    int a; CIN>>a;
/*
    bool flag = false;
    do{
        vector<Shift> shifts = newSolution.getShifts();
        unsigned int shiftIndex, operationIndex;
        shiftIndex = operationIndex = 0;
        flag = false;
        bool breakFlag = false;
        for(int s=0; s<shifts.size(); s++){
            vector<Operation> operations = shifts[s].getOperations();
            operationIndex = 0;
            COUT<<s<<" :"<<" ";
            for(int o=0; o<operations.size(); o++){
                COUT<<operations[o].getPoint()<<" ";
                if(operations[o].getPoint() == customer){
                    shiftIndex = s;
                    operationIndex = o;
                    flag = true;
                    breakFlag = true;
                    break;
                    }
            }
            COUT<<"\n";
            if(breakFlag)
                break;
        }

        if(not flag)
            break;

        if(newSolution.getShifts()[shiftIndex].getOperations().size() <= 1)
            newSolution = this->insertShiftSolution(newSolution, shiftIndex, operationIndex+1, 0, 1.0, 0.0);
        newSolution = this->removeShiftSolution(newSolution, shiftIndex, operationIndex, 1.0, 0.0);
        vector<Shift> newShifts = newSolution.getShifts();

        while(newShifts.size() > shiftIndex+1)
            newShifts.erase(newShifts.end());

        newSolution.setShifts(newShifts);
        newSolution = this->constructSolution(newSolution, 0.1, 0.0, 2, 1.0, 0.0, 0);
        newSolution = this->extendSolution(newSolution, 1.0, 0.0, 0);

        COUT<<"n\n";
        CIN>>a;
    }
    while(flag);
*/

    newSolution = this->constructSolution(newSolution, 0.1, 0.0, 2, 1.0, 0.0, 31, INT_MAX);
    newSolution = this->extendSolution(newSolution, 1.0, 0.0, 31, INT_MAX);

    double oldObjValue = this->computeObjective(newSolution);

    unsigned int shiftIndex, operationIndex;
    shiftIndex = operationIndex = 0;
    newSolution = this->insertShiftSolution(newSolution, 0, 0, 31, 1.0, 0.0);
    double objValue = this->computeObjective(newSolution);
    while(objValue > oldObjValue){

        operationIndex++;
        newSolution = this->insertShiftSolution(newSolution, 0, 0, 31, 1.0, 0.0);


    }

    return newSolution;

}



#endif



