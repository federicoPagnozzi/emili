#ifndef __OPERATION_CPP
#define __OPERATION_CPP

#include <iostream>
#include <sstream>

#include "operation.h"

Operation::Operation()
{
}

unsigned int Operation::getArrival()	{return this->arrival;}
void Operation::setArrival(unsigned int a)	{this->arrival=a;}
unsigned int Operation::getPoint()	{return this->point;}
void Operation::setPoint(unsigned int p)	{this->point=p;}
double Operation::getQuantity()	{return this->quantity;}
void Operation::setQuantity(double q)	{this->quantity=q;}


#endif

