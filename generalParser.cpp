//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include "generalParser.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>

/* modifiers */
#define RO "-ro"
#define IT "-it"
#define RNDSEED "rnds"
#define PRINT_SOLUTION "ps"
#define DEFAULT_TS 10
#define DEFAULT_TI 10
#define DEFAULT_IT 0
#define GIT_COMMIT_NUMBER "aab3f5c3b145b07437c584712d5fa2e74d5320ab"

int tab_level = 0;

void prs::printTab(const char* string)
{
    for(int i=0;i<tab_level; i++)
    {
        std::cout << "  ";
    }

    std::cout << string << std::endl;
}

void prs::incrementTabLevel()
{
    tab_level++;
}

void prs::decrementTabLevel()
{
    tab_level--;
    tab_level= tab_level<0?0:tab_level;
}

void prs::emili_header()
{
    std::cout << "\t ______ __  __ _____ _      _____ " << std::endl;
    std::cout << "\t|  ____|  \\/  |_   _| |    |_   _|" << std::endl;
    std::cout << "\t| |__  | \\  / | | | | |      | |  " << std::endl;
    std::cout << "\t|  __| | |\\/| | | | | |      | |  " << std::endl;
    std::cout << "\t| |____| |  | |_| |_| |____ _| |_ " << std::endl;
    std::cout << "\t|______|_|  |_|_____|______|_____|" << std::endl;
    std::cout << std::endl;    
    std::cout << "commit : " << GIT_COMMIT_NUMBER << std::endl;
}

void prs::info()
{
    std::cout << "Usage:" << std::endl;
    std::cout << std::endl;
    std::cout << "EMILI INSTANCE_FILE_PATH PROBLEM <PROBLEM OPTIONS> [-it|-ro time] [rnds seed] [ps]" << std::endl;
    std::cout << std::endl;
}


void prs::check(char* t,const char* message)
{
    if(t==nullptr)
    {
        prs::info();
        std::cerr <<"PARSING ERROR "<< message << std::endl;
        exit(-1);
    }
}

char* prs::TokenManager::nextToken()
{
    if(currentToken < numberOfTokens)
    {
        char* token = tokens[currentToken];
        currentToken++;
        return token;
    }
    else
    {
        return nullptr;
    }
}

char* prs::TokenManager::peek()
{
    if(currentToken < numberOfTokens)
    {
        char* token = tokens[currentToken];
     //   currentToken++;
        return token;
    }
    else
    {
        return " ";
    }
}

char* prs::TokenManager::operator *()
{
    return peek();
}

void prs::TokenManager::next()
{
    if(currentToken<numberOfTokens)
    {
        currentToken++;
    }
}

void prs::TokenManager::operator ++(int k)
{
    next();
}

char* prs::TokenManager::tokenAt(int i)
{
    if(i < numberOfTokens)
    {
        return tokens[i];
    }
    return nullptr;
}

bool prs::TokenManager::checkToken(std::string &token)
{
    if(currentToken < numberOfTokens)
    {
        if(strcmp(tokens[currentToken],token.c_str())==0)
        {
            currentToken++;
            return true;
        }
    }
        return false;
}

bool prs::TokenManager::checkToken(const char* token)
{
    if(currentToken < numberOfTokens)
    {
        if(strcmp(tokens[currentToken],token)==0)
        {
            currentToken++;
            return true;
        }
    }
        return false;
}

int prs::TokenManager::getInteger()
{
    char* t = peek();
    check(t,"A NUMBER WAS EXPECTED!");
    std::string num(t);
    if( !std::all_of(num.begin(),num.end(),::isdigit))
       check(nullptr,"A NUMBER WAS EXPECTED!");

    int k = atoi(t);
   // std::cout << k << "\n\t";
    next();
    return k;
}

float prs::TokenManager::getDecimal()
{
    char* t = peek();
    check(t,"A DECIMAL NUMBER WAS EXPECTED!");
    std::string num(t);

    if(std::none_of(num.begin(),num.end(),::isdigit))
        check(nullptr,"A DECIMAL NUMBER WAS EXPECTED!");

    float k = atof(t);
    //std::cout << k << "\n\t";
    next();
    return k;
}

bool prs::AlgoBuilder::operator ==(const AlgoBuilder& b)
{
    if(this->availableProblems()==b.availableProblems())
    {
        return true;
    }
    return false;
}


int getTime(prs::TokenManager& tm,int problemSize)
{
        if(tm.checkToken(IT))
        {
            int n = tm.getInteger();
            std::ostringstream oss;
            oss << "Run time secs : " << n;
            //printTab(oss.str().c_str());
            std::cout << oss.str() << std::endl;
            return n;
        }
        else if(tm.checkToken(RO))
        {
            float d = tm.getDecimal();
            float time = d*problemSize;
            int n = floorf(time);
            std::ostringstream oss;
            oss << "Rho = "<< d << " Run time secs : " << n;
            //printTab(oss.str().c_str());
            std::cout << oss.str() << std::endl;
            return n;
        }


    return DEFAULT_IT;
}

int getSeed(prs::TokenManager& tm)
{
    int rnds = 0;
    if(tm.checkToken(RNDSEED))
    {
        rnds = tm.getInteger();
    }
    std::cout << "Random seed : " << rnds << std::endl;
    return rnds;
}


emili::LocalSearch* prs::GeneralParser::parseParams()
{
    tm++;
    tm++;
    char* p = tm.peek();
    if(p != nullptr)
    {
    std::string prob(p);
    for(std::vector< AlgoBuilder*> ::iterator iter= builders.begin(); iter!=builders.end(); ++iter)
    {
        AlgoBuilder* bld = *iter;
        if(bld->isParsable(prob))
        {
            emili::LocalSearch* ls = bld->buildAlgo(tm);
            ls->setSearchTime(getTime(tm,ls->getInitialSolution().getProblem().problemSize()));
            emili::initializeRandom(getSeed(tm));
            if(tm.checkToken(PRINT_SOLUTION))
            {
                emili::set_print(true);
                if(ls->getSearchTime() > 0){
                    std::cout << "The solution will be print at the end of the execution" << std::endl;
                }
            }
            else
            {
                emili::set_print(false);
            }
            return ls;
        }
    }
    std::cout << "------" << std::endl;
    std::cout << "No parser for "<< p << " available" << std::endl;
   }
    std::cerr << "A problem was expected!" << std::endl;
    // info
    prs::info();
 return nullptr;
}

void prs::GeneralParser::registerBuilder(AlgoBuilder* builder)
{
    builders.push_back(builder);
}

void prs::GeneralParser::removeBuilder(AlgoBuilder* builder)
{
    int ind = 0;
    for(;ind<builders.size();ind++)
    {
        if(builder == builders[ind])
        {
            builders.erase(builders.begin()+ind);
            return;
        }
    }
}

