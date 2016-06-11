#ifndef __SOLUTION_CPP
#define __SOLUTION_CPP

#include <iostream>
#include <sstream>

#include "tinyxml.h"
#include "tinystr.h"
#include "instance.h"
#include "solution.h"
#include "operation.h"
#include "instance.h"
#include "driver.h"
#include "trailer.h"
#include "location.h"

#define EPSILON 0.000001



Shift::Shift()
{
}

unsigned int Shift::getIndex()	{return this->index;}
void Shift::setIndex(unsigned int i)	{this->index=i;}
unsigned Shift::getDriver()	{return this->driver;}
void Shift::setDriver(unsigned d)	{this->driver=d;}
unsigned int Shift::getTrailer()	{return this->trailer;}
void Shift::setTrailer(unsigned int t)	{this->trailer=t;}
unsigned int Shift::getStart()	{return this->start;}
void Shift::setStart(unsigned int s)	{this->start=s;}
vector<Operation> Shift::getOperations()	{return this->operations;}
void Shift::setOperations(vector<Operation> o)	{this->operations=o;}
void Shift::setOperation(Operation o, int index)	{this->operations[index]=o;}

double Shift::getShiftQuantity()   {return this->shiftQuantity;}
void Shift::setShiftQuantity(double sq) {this->shiftQuantity = sq;}
double Shift::getShiftDistance()    {return this->shiftDistance;}
void Shift::setShiftDistance(double sd) {this->shiftDistance = sd;}
double Shift::getShiftTime()    {return this->shiftTime;}
void Shift::setShiftTime(double st) {this->shiftTime = st;}

double Shift::getShiftCost()    {return this->shiftCost;}
void Shift::setShiftCost(double sc) {this->shiftCost= sc;}


vector<unsigned int> irpSolution::getRepresentation()  {return this->representation;}
void irpSolution::setRepresentation(vector<unsigned int> representation) {this->representation = representation;}


irpSolution::irpSolution()
{
}

vector<Shift> irpSolution::getShifts()		{return this->shifts;}
void irpSolution::setShifts(vector<Shift> s)	{this->shifts=s;}
void irpSolution::setShift(Shift s, int index)	{this->shifts[index]=s;}
void irpSolution::setPoint(unsigned int point,int sh, int op){
    Shift shift = this->shifts[sh];
    vector<Operation> operations = shift.getOperations();
    Operation operation = operations[op];
    operation.setPoint(point);
    operations[op] = operation;
    shift.setOperations(operations);
    this->setShift(shift, sh);
}

double irpSolution::getTotalQuantity()   {return this->totalQuantity;}
void irpSolution::setTotalQuantity(double tq) {this->totalQuantity = tq;}
double irpSolution::getCost()    {return this->cost;}
void irpSolution::setCost(double c) {this->cost = c;}

vector< vector<double> > irpSolution::getTrailerQuantities()    {return this->trailerQuantities;}
void irpSolution::setTrailerQuantities(vector< vector<double> > tq) {this->trailerQuantities = tq;}
vector< vector<double> > irpSolution::getTankQuantities()   {return this->tankQuantities;}
void irpSolution::setTankQuantities(vector< vector<double> > tq)    {this->tankQuantities = tq;}

void irpSolution::loadSolution(const char *pFilename){
  
      unsigned int index;
      double db;
      TiXmlDocument doc(pFilename);

      if (!doc.LoadFile()) return;
      
      TiXmlHandle hDoc(&doc);
      TiXmlElement* pInnerElem;
      TiXmlHandle hRoot(0);


      // save this for later
      hRoot=TiXmlHandle(pInnerElem);
      
      pInnerElem=hDoc.FirstChild( "IRP_Roadef_Challenge_Output" ).FirstChild("Shifts").FirstChild().Element();
      
      while(pInnerElem!= NULL and strcmp(pInnerElem->Value(), "IRP_Roadef_Challenge_Shift_")==0)
      {
	Shift s;
	TiXmlElement *pShiftElem = pInnerElem->FirstChild()->ToElement();
	
	//Read trailer's index	
	istringstream(pShiftElem->GetText()) >> index;
	s.setIndex(index);	
//	cout<<pShiftElem->Value()<< " " << s.getIndex()<<"\n";
	
	pShiftElem=pShiftElem->NextSiblingElement();
	istringstream(pShiftElem->GetText()) >> index;
	s.setDriver(index);	
//	cout<<pShiftElem->Value()<< " " << s.getDriver()<<"\n";
	
	pShiftElem=pShiftElem->NextSiblingElement();
	istringstream(pShiftElem->GetText()) >> index;
	s.setTrailer(index);	
//	cout<<pShiftElem->Value()<< " " << s.getTrailer()<<"\n";
	
	pShiftElem=pShiftElem->NextSiblingElement();	
	istringstream(pShiftElem->GetText()) >> index;
	s.setStart(index);	
//	cout<<pShiftElem->Value()<< " " << s.getStart()<<"\n";
	
	pShiftElem=pShiftElem->NextSiblingElement();
	TiXmlElement* pOperationElem = pShiftElem->FirstChild()->ToElement();
//	cout<<"        "<<pOperationElem->Value()<<"\n";
	vector<Operation> operations;
	while(pOperationElem!= NULL and strcmp(pOperationElem->Value(), "IRP_Roadef_Challenge_Operation_")==0)
	{
	  Operation o;
	  

	  istringstream(pOperationElem->FirstChild("point")->ToElement()->GetText()) >> index;
	  o.setPoint(index);
//	  cout<<"        "<<pOperationElem->FirstChild("point")->ToElement()->Value()<< " " << o.getPoint()<<"\n";
	  
	  
	  istringstream(pOperationElem->FirstChild("arrival")->ToElement()->GetText()) >> index;
	  o.setArrival(index);	
//	  cout<<"        "<<pOperationElem->FirstChild("arrival")->ToElement()->Value()<< " " << o.getArrival()<<"\n";
	  
	  
	  istringstream(pOperationElem->FirstChild("Quantity")->ToElement()->GetText()) >> db;
	  o.setQuantity(db);	
//	  cout<<"        "<<pOperationElem->FirstChild("Quantity")->ToElement()->Value()<< " " << o.getQuantity()<<"\n";
	  
	  pOperationElem=pOperationElem->NextSiblingElement();
	  
	  operations.push_back(o);
	}
	
	s.setOperations(operations);
	this->shifts.push_back(s);
	
	pInnerElem=pInnerElem->NextSiblingElement();
      }
      
}



void irpSolution::saveSolution(string pFilename){
  
  ofstream file;
  file.open (pFilename);
  file.precision(15);
  file << "<IRP_Roadef_Challenge_Output>\n<Shifts>";
  
  for(int s=0; s<this->shifts.size(); s++){
    Shift shift = this->shifts[s];
    file << "<IRP_Roadef_Challenge_Shift_>";
    file <<"\n<index>"<<shift.getIndex()<<"</index>";
    file <<"\n<driver>"<<shift.getDriver()<<"</driver>";
    file <<"\n<trailer>"<<shift.getTrailer()<<"</trailer>";
    file <<"\n<start>"<<shift.getStart()<<"</start>";
    
    file <<"\n<operations>";
    vector<Operation> operations = shift.getOperations();
    for(int o=0; o<operations.size(); o++){
      Operation operation = operations[o];
      file <<"\n<IRP_Roadef_Challenge_Operation_>";
      file <<"\n<point>"<<operation.getPoint()<<"</point>";
      file <<"\n<arrival>"<<operation.getArrival()<<"</arrival>";
      file <<"\n<Quantity>"<<operation.getQuantity()<<"</Quantity>";
      file <<"\n</IRP_Roadef_Challenge_Operation_>";
    }
    file <<"\n</operations>";
      
    file << "\n</IRP_Roadef_Challenge_Shift_>";
  }
  
  
  file << "\n</Shifts>\n</IRP_Roadef_Challenge_Output>";
  file.close();
}

void irpSolution::fromSolutionToRepresentation(irpSolution solution){

    vector<unsigned int> representation;
    for(int s=0; s<solution.getShifts().size(); s++)
        for(int o=0; o<solution.getShifts()[s].getOperations().size(); o++)
            if(solution.getShifts()[s].getOperations()[o].getPoint() != 1)
                representation.push_back(solution.getShifts()[s].getOperations()[o].getPoint());
    this->representation = representation;
}

/*
void irpSolution::readSolutionFromInstance(Instance instance){
  
  vector<Shift>::iterator it;
  
  cout<<"here "<<this->shifts.at(0).getOperations().at(0).getPoint().getIndex()<<"\n";
  for( it = this->shifts.begin(); it < this->shifts.end(); ++it) { 
    cout<<"Shift number "<<it->getIndex()<<"\n";
    it->setDriver(instance.getDrivers().at(it->getDriver().getIndex()));
//    cout<<it->getDriver().getIndex() <<"	"<<it->getDriver().getMaxDrivingDuration()<<"\n";
    it->setTrailer(instance.getTrailers().at(it->getTrailer().getIndex()));
    
     
    for( int i=0; i< it->getOperations().size(); i++) { 
      vector<Operation> o=it->getOperations();
      Customer c;
      vector<Operation> o2=it->getOperations();
      cout<<o.at(i).getPoint().getIndex()<< " ";
      if(o.getPoint().getIndex()==0)
	o2.push_back(instance->getBase());
      else if(o.getPoint().getIndex()==1)
	o2.push_back(instance->getSource());
      else
	o2.push_back(instance.getCustomers().at(o.at(i).getPoint().getIndex()));
    
    }
    cout<<"\n";
    
  }
}
*/
#endif

