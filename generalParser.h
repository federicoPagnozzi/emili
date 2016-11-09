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
/* This method prints the big emili at the beginning of the execution.
 */
void emili_header();
/* info prints the basic informations to use EMILI.
 * It does not print specific informations about the AlgoBuilder,
 * those have to be provided by the Algobuilder class.
 */
void info();
/* Checks that t is not nullptr otherwise
 * it prints message and halts the execution.
 */
void check(char* t,const char* message);
/* printTab is used to print in a more readable fashion
 * messages printed during the parsing process.
 * when incrementTabLevel is called an additional \t is added ad the beginning
 * of the line before string. decrementTabLevel reduces the number of \t
 * displayed.
 */
void printTab(const char* string);
void incrementTabLevel();
void decrementTabLevel();

/* This class implements the Token Manager used
 * to parse the arguments given to EMILI
 */
class TokenManager
{
protected:
    char** tokens;
    int numberOfTokens;
    int currentToken;
    int previousCurrentToken;
public:
    TokenManager(char** tokens,int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(0),previousCurrentToken(0) { }
    /*
     * comsumes a token and return it
        */
    char* nextToken();
    /*
     * returns a token without consuming it
     */
    char* peek();
    // same as peek
    char* operator *();
    /*
     * consumes a token
    */
    void next();
    // same as next
    void operator ++(int);
    /*
     * This method returns an integer
     * If current token can be parsed as an integer.
     * otherwise it will show an error and end the execution.
     */
    int getInteger();
    /*
     * This method returns a float
     * If current token can be parsed as a float.
     * otherwise it will show an error and end the execution.
     */
    float getDecimal();
    /*
     * check if the parameter matches the current token
     *  if it does the token is consumed and it returns true.
    */
    bool checkToken(std::string& token);
    bool checkToken(const char* );
    /*
     * returns the token at position i
    */
    char* tokenAt(int i);
    /*
     * return the position of the string in the token list if present, -1 otherwise
    */
    int seek(const char*);
    /*
     * moves the current token to position i and return true.
     * returns false if the position is less than zero or more than the number of tokens
     *
    */
    bool move(int i);
    /*
     * restores the current token index to before the last move operation
    */
    void restore();
    /*
    * returns true if there are more tokens to parse
    */
    bool hasMoreTokens();

};

/* AlgoBuilder instantiate an algorithm starting from the parameters given at running time.
 *
 */
class AlgoBuilder
{
protected:
    /*This method should return a string in which are listed all the problems supported */
    virtual std::string availableProblems() const{ return std::string("Iamabstract!");}
public:
    /*This method should return true if the object is capable of building an algorithm
      to solve problem */
    virtual bool isParsable(std::string& problem)=0 ;
    /* This method should return an algorithm ready to be run.
     * The algorithm is built by parsing the configuration provided a run time
     * and contained in tm.
     */
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm) {return nullptr;}
    /* this method should return a description of all the differents components
     * that are supported by the object.
     */
    virtual std::string info() {return std::string("Iamabstract!");}
    /* Overloading of the == operator for Algobuilder.
     * By default two Algobuilder are considered equal if they support the same problem(s)
     */
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
