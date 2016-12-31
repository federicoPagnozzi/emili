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
#define DEFAULT_TS 10
#define DEFAULT_TI 10
#define DEFAULT_IT 0



int tab_level = 0;

void prs::printTab(const std::string& string) {
    printTab(string.c_str());
}

void prs::printTab(const char* string)
{
    for(int i=0;i<tab_level; i++)
        std::cout << "| ";

    std::cout << string << std::endl;
}

void prs::incrementTabLevel()
{
    tab_level++;
}

void prs::decrementTabLevel()
{
    tab_level--;
    if(tab_level < 0)
        tab_level = 0;
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
}

void prs::info()
{
    std::cout << "Usage:" << std::endl;
    std::cout << std::endl;
    std::cout << "EMILI INSTANCE_FILE_PATH PROBLEM <PROBLEM OPTIONS> [-it|-ro time] [rnds seed]" << std::endl;
    std::cout << std::endl;
}


void prs::check(const char* t, const char* message)
{
    if(t == nullptr)
    {
        // prs::info();
        std::cerr << "PARSING ERROR : " << message << std::endl;
        exit(-1);
    }
}

const char* prs::TokenManager::nextToken()
{
    if(currentToken < numberOfTokens)
    {
        const char* token = tokens[currentToken];
        currentToken++;
        return token;
    }
    else
    {
        return nullptr;
    }
}

/**
 * @brief prs::TokenManager::peek
 * @return current token or " "
 */
const char* prs::TokenManager::peek() const
{
    if(currentToken < numberOfTokens)
        return tokens[currentToken];
    else
        return " ";
}

const char* prs::TokenManager::operator *() const
{
    return peek();
}

void prs::TokenManager::next()
{
    if(currentToken < numberOfTokens)
        currentToken++;
}

bool prs::TokenManager::hasNext() const {
    return currentToken < numberOfTokens;
}

void prs::TokenManager::operator ++(int k)
{
    next();
}

const char* prs::TokenManager::tokenAt(int i) const
{
    if(i < numberOfTokens)
        return tokens[i];

    return nullptr;
}

bool prs::TokenManager::checkToken(std::string const& token)
{
    if(peekIs(token)) {
        currentToken++;
        return true;
    } else {
        return false;
    }
}

bool prs::TokenManager::checkToken(const char* token)
{
    if(peekIs(token)) {
        currentToken++;
        return true;
    } else {
        return false;
    }
}

bool prs::TokenManager::peekIs(const char* token) const {
    return currentToken < numberOfTokens && strcmp(tokens[currentToken], token) == 0;
}

bool prs::TokenManager::peekIs(std::string const& token) const {
    return currentToken < numberOfTokens && tokens[currentToken] == token;
}

bool prs::TokenManager::checkInteger(int & res) {
    if(currentToken < numberOfTokens) {
        std::istringstream iss(tokens[currentToken]);
        iss >> res;
        if(iss) {
            currentToken++;
            return true;
        }
    }
    return false;
}

bool prs::TokenManager::checkDecimal(float & res) {
    if(currentToken < numberOfTokens) {
        std::istringstream iss(tokens[currentToken]);
        iss >> res;
        if(iss) {
            currentToken++;
            return true;
        }
    }
    return false;
}

bool prs::TokenManager::checkDecimal(double & res) {
    if(currentToken < numberOfTokens) {
        std::istringstream iss(tokens[currentToken]);
        iss >> res;
        if(iss) {
            currentToken++;
            return true;
        }
    }
    return false;
}

namespace {
float toDecimal(std::string const& s) throw(prs::ParsingError) {
    float k;
    std::istringstream iss(s);
    iss >> k;
    if(!iss)
        throw prs::ParsingError(std::string("Decimal expected, '") + s + "' found");
    return k;
}
}

int prs::TokenManager::getInteger() throw(ParsingError)
{
    int k;
    std::istringstream iss(peek());
    iss >> k;
    if(!iss)
        throw ParsingError(std::string("Int expected, '") + peek() + "' found");

    next();
    return k;
}

float prs::TokenManager::getDecimal() throw(ParsingError)
{
    float k = toDecimal(peek());
    next();
    return k;
}

float prs::TokenManager::getPercent(float N) throw(ParsingError)
{
    if(checkToken("percent") || checkToken("%")) {
        return getDecimal() * N / 100;
    } else {
        try {
            return getDecimal();
        } catch(ParsingError&) {
            std::string p = peek(), tok;
            if(p.substr(p.size() - 1) == "%") {
                tok = p.substr(0, p.size() - 1);
            } else if(p.substr(0, 1) == "%") {
                tok = p.substr(1);
            } else {
                tok = peek();
            }
            float k = toDecimal(tok);
            next();
            return k * N / 100;
        }
    }
}

bool prs::AlgoBuilder::operator ==(const AlgoBuilder& b)
{
    return this->availableProblems() == b.availableProblems();
}

emili::LocalSearch* prs::GeneralParser::parseParams()
{
    tm++; // tokenAt(0) == argv[0]
    tm++; // tokenAt(1) == instancePath
    std::string problem = tm.peek(); // tokenAt(2) == problem
    for(AlgoBuilder* builder : builders) {
        if(builder->isParsable(problem))
        {
            emili::LocalSearch* ls = builder->buildAlgo(tm);
            emili::Problem* instance = builder->getInstance();
            if(! instance)
                instance = &ls->getInitialSolution().getProblem();

            noSearch = false;
            int it = DEFAULT_IT;
            int seed = 0;
            for(;;) {
                if(tm.checkToken("no-search"))
                    noSearch = true;
                else if(tm.checkToken(IT))
                    it = tm.getInteger();
                else if(tm.checkToken(RO))
                    it = floorf(tm.getDecimal() * instance->problemSize());
                else if(tm.checkToken(RNDSEED))
                    seed = tm.getInteger();
                else
                    break;
            }
            if(seed == 0)
                seed = emili::getRandomSeedFromRandom();
            std::cout << "Run time secs : " << it << std::endl;
            std::cout << "Random seed : " << seed << std::endl;
            ls->setSearchTime(it);
            emili::initializeRandom(seed);
            return ls;
        }
    }
    std::vector<std::string> avail;
    for(AlgoBuilder* bld : builders) {
        auto L = bld->availableProblemsList();
        avail.insert(avail.end(), L.begin(), L.end());
    }
    throw ErrorExpected(tm, "problem", avail);
}

void prs::GeneralParser::registerBuilder(AlgoBuilder* builder)
{
    builders.push_back(builder);
}

void prs::GeneralParser::removeBuilder(AlgoBuilder* builder)
{
    auto it = std::find(builders.begin(), builders.end(), builder);
    if(it != builders.end())
        builders.erase(it);
}

