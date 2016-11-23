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
/** This method prints the big emili at the beginning of the execution.
 */
void emili_header();
/** info prints the basic informations to use EMILI.
 * It does not print specific informations about the AlgoBuilder,
 * those have to be provided by the Algobuilder class.
 */
void info();
/** Checks that t is not nullptr otherwise
 * it prints message and halts the execution.
 */
void check(char* t,const char* message);
/** printTab is used to print in a more readable fashion
 * messages printed during the parsing process.
 * when incrementTabLevel is called an additional \ t is added ad the beginning
 * of the line before string. decrementTabLevel reduces the number of \ t
 * displayed.
 */
void printTab(const char* string);
void incrementTabLevel();
void decrementTabLevel();

/** This class implements the Token Manager used
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
    /**
     * comsumes a token and return it
        */
    char* nextToken();
    /**
     * returns a token without consuming it
     */
    char* peek();
    /** same as peek */
    char* operator *();
    /**
     * consumes a token
    */
    void next();
    /** same as next*/
    void operator ++(int);
    /**
     * This method returns an integer
     * If current token can be parsed as an integer.
     * otherwise it will show an error and end the execution.
     */
    int getInteger();
    /**
     * This method returns a float
     * If current token can be parsed as a float.
     * otherwise it will show an error and end the execution.
     */
    float getDecimal();
    /**
     * check if the parameter matches the current token
     *  if it does the token is consumed and it returns true.
    */
    bool checkToken(std::string& token);
    bool checkToken(const char* );
    /**
     * returns the token at position i
    */
    char* tokenAt(int i);
    /**
     * return the position of the string in the token list if present, -1 otherwise
    */
    int seek(const char*);
    /**
     * moves the current token to position i and return true.
     * returns false if the position is less than zero or more than the number of tokens
     *
    */
    bool move(int i);
    /**
     * restores the current token index to before the last move operation
    */
    void restore();
    /**
    * returns true if there are more tokens to parse
    */
    bool hasMoreTokens();

};

/** AlgoBuilder instantiate an algorithm starting from the parameters given at running time.
 * This class ( and also GeneralParser ) is part of the old system for parsing the command line
 * arguments. If you are implementing components for a new problem or planning to add other to
 * an already imlpemented one you shoudl check the new system
 * !these classes will be removed soon from the main branch!
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
    /** This method should return an algorithm ready to be run.
     * The algorithm is built by parsing the configuration provided a run time
     * and contained in tm.
     */
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm) {return nullptr;}
    /** this method should return a description of all the differents components
     * that are supported by the object.
     */
    virtual std::string info() {return std::string("Iamabstract!");}
    /** Overloading of the == operator for Algobuilder.
     * By default two Algobuilder are considered equal if they support the same problem(s)
     */
    virtual bool operator ==(const AlgoBuilder& b);
};

/** This class as well as Algobuilder are part of the old system for parsing the command line
* arguments. If you are implementing components for a new problem or planning to add other to
* an already imlpemented one you shoudl check the new system
* !these classes will be removed soon from the main branch!*/
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

/**
 * TYPES TABLE
 * In this "table" are defined the type values for the
 * types defined in emilibase.h
 */
/**    ---------------------------------------------------
 *    | Component                            | Value
      ---------------------------------------------------*/
#define COMPONENT_EMPTY                         0xEE //This type model the empty symbol
#define COMPONENT_NULL                          0xFF //Something was expected but nothing was found
#define COMPONENT_ALGORITHM                     0xA0
#define COMPONENT_INITIAL_SOLUTION_GENERATOR    0xB1
#define COMPONENT_TERMINATION_CRITERION         0xB2
#define COMPONENT_NEIGHBORHOOD                  0xB3
#define COMPONENT_TABU_TENURE                   0xB4
#define COMPONENT_NEIGHBORHOOD_OR_EMPTY         0xBE
#define COMPONENT_PERTURBATION                  0xC1
#define COMPONENT_ACCEPTANCE                    0xC2


/**
 * A Component rapresent one of the parts in which the algorithm is divided.
 * It can be of many types ( see types table ) from an Algorithm
 * to a Termination criterion
 */
class Component
{
protected:
    /** The type of the component (see types table) */
    int type;
    /** The pointer to the actual component*/
    void* rawComponent;
    /** The string token that represents this component (not used)*/
    char* token;
public:
    Component(int type,void* rawData):type(type),rawComponent(rawComponent) { }
    /** The default component is NULL*/
    Component():type(COMPONENT_NULL),rawComponent(nullptr) { }
    virtual Component& operator=(const Component& a);
    /** Return true if the component is of type type*/
    bool is(int type) {return this->type == type;}
    int getType() {return type;}
    void setType(int type) {this->type = type;}
    void setRawComponent(void* rawc) {this->rawComponent = rawc;}
    /** Casts rawCompoent to type T* */
    template<typename T>
    T* get(){ return (T*) rawComponent;}
};

class GeneralParserE;
class Builder
{
protected:
    /**/
    GeneralParserE& gp;
    TokenManager& tm;
    virtual Component retrieveComponent(int type,bool allowEmpty=false);
    virtual std::vector<emili::Neighborhood*> buildNeighborhoodVector();
public:
    Builder(GeneralParserE& generalParser,TokenManager& tokenManager):gp(generalParser),tm(tokenManager) { }
    virtual emili::Problem* openInstance() {return nullptr;}
    virtual bool isCompatibleWith(char* problem_definition)=0;
    virtual bool isCompatibleWith(std::string& problem_definition) {return isCompatibleWith((char*)problem_definition.c_str());}
    virtual bool canOpenInstance(char* problem_definition) {return false;}
    virtual bool canOpenInstance(std::string& problem_definition) {return canOpenInstance((char*)problem_definition.c_str());}
    virtual std::string typeName(int type){return std::string();}
    virtual Component buildComponent(int type);
    virtual emili::LocalSearch* buildAlgo() {return nullptr;}
    virtual emili::InitialSolution* buildInitialSolution(){return nullptr;}
    virtual emili::Neighborhood* buildNeighborhood(){return nullptr;}
    virtual emili::Termination* buildTermination(){return nullptr;}
    virtual emili::Perturbation* buildPerturbation(){return nullptr;}
    virtual emili::Acceptance* buildAcceptance(){return nullptr;}
    virtual emili::TabuMemory* buildTabuTenure(){return nullptr;}
    virtual std::string cmd_info(){return std::string("");}
};

class GeneralParserE: public GeneralParser
{
protected:
    std::vector< Builder* > allbuilders;
    std::vector< Builder* > activeBuilders;
    emili::Problem* instance;
    virtual std::string typeName(int type);
public:
    GeneralParserE(char** tokens,int numberOfTokens):GeneralParser(tokens,numberOfTokens) { }
    GeneralParserE(TokenManager tokenmanager):GeneralParser(tokenmanager) { }
    virtual emili::LocalSearch* parseParams();
    virtual void addBuilder(Builder* b);
    virtual emili::Problem* getInstance() {return instance;}
    virtual void setInstance(emili::Problem* ins) {instance=ins;}
    virtual TokenManager& getTokenManager() {return tm;}
    virtual Component buildComponent(int type);
    virtual void fatalError(int received_type,int expected_type);    
};

class EmBaseBuilder: public Builder
{
public:
    EmBaseBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char *problem_definition) { return true;}
    virtual emili::LocalSearch* buildAlgo();
    virtual emili::Termination* buildTermination();
    virtual emili::Perturbation* buildPerturbation();
    virtual emili::Acceptance* buildAcceptance();
};

}

#endif
