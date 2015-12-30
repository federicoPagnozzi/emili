#ifndef __LOCATION_CPP
#define __LOCATION_CPP

#include "location.h"


unsigned int Customer::getIndex()	{return this->index;}
void Customer::setIndex(unsigned int ind)	{this->index=ind;}

unsigned int Customer::getSetupTime()	{return this->setupTime;}
void Customer::setSetupTime(unsigned int st)	{this->setupTime=st;}

Customer::Customer()	{}
vector<unsigned int> Customer::getAllowedTrailers()	{return this->allowedTrailers;}
void Customer::setAllowedTrailers(vector<unsigned int> at)	{this->allowedTrailers=at;}
vector<double> Customer::getForecast()	{return this->forecast;}
void Customer::setForecast(vector<double> f)	{this->forecast=f;}
double Customer::getCapacity()	{return this->capacity;}
void Customer::setCapacity(double c)	{this->capacity=c;}
double Customer::getInitialTankQuantity()	{return this->initialTankQuantity;}
void Customer::setInitialTankQuantity(double itq)	{this->initialTankQuantity=itq;}
double Customer::getSafetyLevel()	{return this->safetyLevel;}
void Customer::setSafetyLevel(double sl)	{this->safetyLevel=sl;}

#endif


