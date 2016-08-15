//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#ifndef GENERALPARSER_H
#define GENERALPARSER_H
#include "emilibase.h"

#include <exception>
#include <string>
#include <vector>
#include <sstream>

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

struct ParsingError : std::exception {
    std::string msg;
    ParsingError() {}
    ParsingError(std::string msg_) : msg(msg_) {}

    const char* what() const throw() override {
        return msg.c_str();
    }
};

class TokenManager
{
protected:
    const char** tokens;
    int numberOfTokens;
    int currentToken;
public:
    TokenManager(const char** tokens, int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(0) { }

    bool hasNext();
    void next();
    void operator ++(int); // same as next
    const char* nextToken(); // skip one token and return it, or null if no token left

    const char* peek(); // return current token or " "
    const char* operator *(); // same as peek

    int getInteger() throw(ParsingError);
    float getDecimal() throw(ParsingError);

    bool checkToken(const std::string &token);
    bool checkToken(const char*);

    const char* tokenAt(int i);
    bool peekIs(const char*);
    bool peekIs(const std::string &);

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

struct ErrorExpected : ParsingError {
    ErrorExpected(prs::TokenManager& tm, std::string name, std::vector<std::string> const& tokens) {
        std::ostringstream oss;
        oss << "'" << tm.peek() << "' -> ERROR a " << name << " is expected : ";
        if(tokens.size() == 0) {
            oss << "<>";
        } else {
            auto it = tokens.begin();
            oss << "<" << *it++;
            while(it != tokens.end())
                oss << " | " << *it++;
            oss << ">";
        }
        msg = oss.str();
    }
};

class AlgoBuilder
{
public:
    virtual std::string availableProblems() const{ return std::string("Iamabstract!");}
    virtual std::vector<std::string> availableProblemsList() const { return {}; }
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
    const char* what() const throw() {
        return "NoSearch";
    }
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
    bool noSearch = false;
};
}

#endif
