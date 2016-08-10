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
     * comsumes the token and return it
        */
    char* nextToken();
    /*
     * returns the token without consuming it
     */
    char* peek();
    // same as peek
    char* operator *();
    /*
     * consumes the token
    */
    void next();
    // same as next
    void operator ++(int);
    int getInteger();
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

class AlgoBuilder
{
protected:
    virtual std::string availableProblems() const{ return std::string("Iamabstract!");}
public:
    virtual bool isParsable(std::string& problem)=0 ;
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm) {return nullptr;}
    virtual std::string info() {return std::string("Iamabstract!");}
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
