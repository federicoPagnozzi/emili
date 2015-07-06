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
public:
    TokenManager(char** tokens,int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(0) { }
    char* nextToken();
    char* peek();
    // same as peek
    char* operator *();
    void next();
    // same as next
    void operator ++(int);
    int getInteger();
    float getDecimal();
    bool checkToken(std::string& token);
    bool checkToken(const char* );
    char* tokenAt(int i);

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
