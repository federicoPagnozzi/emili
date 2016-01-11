//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#ifndef GENERALPARSER_H
#define GENERALPARSER_H
#include "emilibase.h"

namespace prs
{
void emili_header();
void info();
void check(char* t,const char* message);
void printTab(const char* string);
void incrementTabLevel();
void decrementTabLevel();

/*
 * This class implements a TokenManager for the parsing of the arguments
   It does all the dirty work for number parsing and string comparing.
*/
class TokenManager
{
protected:
    char** tokens;
    int numberOfTokens;
    int currentToken;
public:
    TokenManager(char** tokens,int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(0) { }
    //takes the next token from the queue and returns it
    char* nextToken();
    //returns the next token without removing it from the queue
    char* peek();
    // same as peek
    char* operator *();

    //remove the next token from the queue
    void next();
    // same as next
    void operator ++(int);

    /* remove the next token from the queue and returns it as an integer value.
     If the conversion to integer fails the method prints on standard output an error message and
     stops the execution.*/
    int getInteger();
    /* remove the next token from the queue and returns it as a floating point value.
     If the conversion fails the method prints on standard output an error message and
     stops the execution.*/
    float getDecimal();
    /*Checks that the given string is equal to the next token.
      If there is a match the token is removed from the queue and the methods returns true.
      Otherwise the method returns false and it does not remove the token from the queueu.
    */
    bool checkToken(std::string& token);
    /*Same as checkToken with a std::string*/
    bool checkToken(const char* );
    /*returns the token at the index i*/
    char* tokenAt(int i);

};

class AlgoBuilder
{
protected:
    /*Should return all the different problem for which this class can build an algorithm*/
    virtual std::string availableProblems() const{ return std::string("Iamabstract!");}
public:
    /*This method tells the general parser if this class can parse the arguments for problem*/
    virtual bool isParsable(std::string& problem)=0 ;
    /*
        This method is called by the general parser to build the algorithm.
        It must use the TokenManager to parse the arguments and build the algorithm accordigly.
    */
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm) {return nullptr;}
    /*This method should return human readable informations to let the user know how to instantiate
      an algorithm for this problem : which kind of metaheuristics are supported (ILS, TABU, SA etc...),
      which kind of neighborhood, termination criteria, acceptance , perturbations and so on...
    */
    virtual std::string info() {return std::string("Iamabstract!");}
    /* */
    virtual bool operator ==(const AlgoBuilder& b);
};

class GeneralParser
{
protected:
    std::vector< AlgoBuilder* > builders;
    TokenManager tm;
public:        
    GeneralParser(char** tokens,int numberOfTokens):tm(tokens,numberOfTokens) { }
    GeneralParser(TokenManager tokenmanager):tm(tokenmanager) { }
    virtual emili::LocalSearch* parseParams();
    virtual void registerBuilder(AlgoBuilder* builder);
    virtual void removeBuilder(AlgoBuilder* builder);
};
}

#endif
