//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#ifndef GENERALPARSER_H
#define GENERALPARSER_H
#include "emilibase.h"

#include <exception>

namespace prs
{
void emili_header();
void info();
void check(const char *t, const char* message);
void printTab(const char* string);
void printTab(const std::string& string);
void incrementTabLevel();
void decrementTabLevel();

struct TabLevel {
    TabLevel() { incrementTabLevel(); }
    ~TabLevel() { decrementTabLevel(); }
};

class TokenManager
{
protected:
    const char** tokens;
    int numberOfTokens;
    int currentToken;
public:
    TokenManager(const char** tokens,int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(0) { }
    const char *nextToken();
    const char *peek();
    // same as peek
    const char *operator *();
    void next();
    // same as next
    void operator ++(int);
    int getInteger();
    float getDecimal();
    bool checkToken(const std::string &token);
    bool checkToken(const char* );
    const char* tokenAt(int i);

    /**
     * @brief if(the next token is an integer) { res = the integer; return true } else { return false }
     * @param res the location to write the integer
     * @return true if an integer was set to res
     */
    bool checkInteger(int &res);

    /**
     * @brief if(the next token is a decimal) { res = the decimal; return true } else { return false }
     * @param res the location to write the decimal
     * @return true if an decimal was set to res
     */
    bool checkDecimal(float &res);

    /**
     * @brief if(the next token is a decimal) { res = the decimal; return true } else { return false }
     * @param res the location to write the decimal
     * @return true if an decimal was set to res
     */
    bool checkDecimal(double &res);
};

class AlgoBuilder
{
protected:
    virtual std::string availableProblems() const{ return std::string("Iamabstract!");}
public:
    virtual bool isParsable(std::string& problem)=0 ; // should be const std::string &
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm) {return nullptr;}
    virtual std::string info() {return std::string("Iamabstract!");}
    virtual bool operator ==(const AlgoBuilder& b);

    /**
     * may return null
     */
    virtual emili::Problem* getInstance() {
        return nullptr;
    }
};

class NoSearch : public std::exception {

};

class GeneralParser
{
protected:
    std::vector< AlgoBuilder* > builders;
    TokenManager tm;
public:        
    GeneralParser(const char** tokens, int numberOfTokens):tm(tokens,numberOfTokens) { }
    GeneralParser(TokenManager tokenmanager):tm(tokenmanager) { }
    virtual emili::LocalSearch* parseParams();
    virtual void registerBuilder(AlgoBuilder* builder);
    virtual void removeBuilder(AlgoBuilder* builder);
};
}

#endif
