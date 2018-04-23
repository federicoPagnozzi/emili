//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#ifndef GENERALPARSER_H
#define GENERALPARSER_H
#include "emilibase.h"
/**
 * All the classes that are involved in the parsing of the command line belongs to this namespace
 */
namespace prs
{
/**
 * @brief emili_header
 * This method prints the big emili at the beginning of the execution.
 */
void emili_header();
/**
 * @brief info
 * info prints the basic informations to use EMILI.
 * It does not print specific informations about the AlgoBuilder,
 * those have to be provided by the Algobuilder class.
 */
void info();
/**
 * @brief check
 *  Checks that t is not nullptr otherwise
 * it prints message and halts the execution.
 * @param t
 * @param message
 */
void check(char* t,const char* message);

/**
 * @brief printTab
 * printTab is used to print in a more readable fashion
 * messages printed during the parsing process.
 * when incrementTabLevel is called an additional \ t is added ad the beginning
 * of the line before string. decrementTabLevel reduces the number of \ t
 * displayed.
 * @param string
 */
void printTab(const char* string);
/**
 * @brief printTab
 * printTabPlusOne is used to print messages one tabLevel more than printTab
 * when incrementTabLevel is called an additional \ t is added ad the beginning
 * of the line before string. decrementTabLevel reduces the number of \ t
 * displayed.
 * @param string
 */
void printTabPlusOne(const char* string);

int getTabLevel();

template <typename T >
void printTabPlusOne(const char* string,T value)
{
    int tab_level = getTabLevel();
    for(int i=0;i<=tab_level; i++)
    {
        std::cout << "  ";
    }
    std::cout << string << " : " << value << std::endl;
}
/**
 * @brief incrementTabLevel
 *          Increments the number of tabs added at the beginning of a line by printTab
 */
void incrementTabLevel();
/**
 * @brief decrementTabLevel
 *          Decrements the number of tabs added at the beginning of a line by printTab
 */
void decrementTabLevel();
/**
 * @brief The TokenManager class
 * This class implements the Token Manager used
 * to parse the arguments given to EMILI
 */
class TokenManager
{
protected:
    /**
     * @brief tokens
     *          The pointer to the command line arguments
     */
    char** tokens;
    /**
     * @brief numberOfTokens
     *          Total number of tokens.
     */
    int numberOfTokens;
    /**
     * @brief currentToken
     *          index to the current token.
     */
    int currentToken;
    /**
     * @brief previousCurrentToken
     *          index to the previuos value assumed by currentToken
     */
    int previousCurrentToken;
    /**
     * @brief empty
     *        returned in case there is no token
     */
    char empty[3] = {'"',' ','"'};
public:
    /**
     * @brief TokenManager
     *          The constructor needs a pointer to the commandline arguments and the number of
     *          elements in the array.
     * @param tokens
     *           Pointer to commandline arguments array.
     * @param numberOfTokens
     *           Length of the array.
     */
    TokenManager(char** tokens,int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(0),previousCurrentToken(0) { }
    /**
     * @brief nextToken
     *           comsumes a token and return it
     * @return
     *          char string representing the token
     */
    char* nextToken();
    /**
     * @brief peek
     *      returns a token without consuming it
     * @return
     *          char string representing the token
     */
    char* peek();
    /**
     * @brief operator *
     *      Works like peek.
     *      returns a token without consuming it
     * @return
     *          char string representing the token
     */
    char* operator *();
    /**
     * @brief next
     * consumes a token
     * @return
     *          char string representing the token
     */
    void next();
    /**
     * @brief operator ++
     * consumes a token
     * @return
     *          char string representing the token
     */
    void operator ++(int);
    /**
     * @brief getInteger
     *
     * If current token can be parsed as an integer.
     * otherwise it will show an error and end the execution.
     * @return
     * This method returns an integer
     */
    int getInteger();
    /**
     * @brief getDecimal
     * If current token can be parsed as a float.
     * otherwise it will show an error and end the execution.
     * @return
     * This method returns a float
     */
    float getDecimal();
    /**
     * @brief checkToken
     *  check if token matches currentToken
     *  if it does the token is consumed and it returns true.
     * @param token
     * @return
     *  true if token matches currentToken, false otherwise
     */
    bool checkToken(std::string& token);
    /**
     * @brief checkToken
     *  check if token matches currentToken
     *  if it does the token is consumed and it returns true.
     * @param token
     * @return
     *  true if token matches currentToken, false otherwise
     */
    bool checkToken(const char* );
    /**
     * @brief tokenAt
     * @param i
     * @return
     *          returns the token at position i
     */
    char* tokenAt(int i);    
    /**
     * @brief seek
     * @return
     * return the position of the string in the token list if present, -1 otherwise
     */
    int seek(const char*);    
    /**
     * @brief move
     * moves the current token to position i and return true.
     * @param i
     * @return
     * returns false if the position is less than zero or more than the number of tokens
     */
    bool move(int i);
    /**
     * @brief restore
     * restores the current token index to before the last move operation
     */
    void restore();
    /**
     * @brief hasMoreTokens
     * @return
     * returns true if there are more tokens to parse
     */
    bool hasMoreTokens();

};

/**
 * @brief The AlgoBuilder class
 * AlgoBuilder instantiate an algorithm starting from the parameters given at running time.
 * This class ( and also GeneralParser ) is part of the old system for parsing the command line
 * arguments. If you are implementing components for a new problem or planning to add other to
 * an already imlpemented one you shoudl check the new system
 * !these classes will be removed soon from the main branch!
 */
class AlgoBuilder
{
protected:
    /**
     * @brief availableProblems
     * This method should return a string in which are listed all the problems supported
     * @return
     */
    virtual std::string availableProblems() const{ return std::string("Iamabstract!");}    
public:
    /**
     * @brief isParsable
     * This method should return true if the object is capable of building an algorithm
     *  to solve problem
     * @param problem
     * @return
     */
    virtual bool isParsable(std::string& problem)=0 ;
    /**
     * @brief buildAlgo
     * This method should return an algorithm ready to be run.
     * The algorithm is built by parsing the configuration provided a run time
     * and contained in tm.
     * @param tm
     * @return
     */
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm) {return nullptr;}
    /**
     * @brief info
     * this method should return a description of all the differents components
     * that are supported by the object.
     */
    virtual std::string info() {return std::string("Iamabstract!");}
    /**
     * @brief operator ==
     * Overloading of the == operator for Algobuilder.
     * By default two Algobuilder are considered equal if they support the same problem(s)
     */
    virtual bool operator ==(const AlgoBuilder& b);
};

/**
 * @brief The GeneralParser class
 * This class as well as Algobuilder are part of the old system for parsing the command line
 * arguments. If you are implementing components for a new problem or planning to add other to
 * an already imlpemented one you shoudl check the new system
 * !these classes will be removed soon from the main branch!
 */
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
#define COMPONENT_SHAKE                         0xC3
#define COMPONENT_NEIGHBORHOOD_CHANGE           0xC4
#define COMPONENT_PROBLEM                       0x99 //Request to load a problem


/**
 * @brief Component class
 * A Component rapresent one of the parts in which the algorithm is divided.
 * It can be of many types ( see types table ) from an Algorithm
 * to a Termination criterion
 */
class Component
{
protected:    
    /**
     * @brief type
     * The type of the component (see types table)
      ---------------------------------------------------
      | Component                            | Value
      ---------------------------------------------------
       COMPONENT_EMPTY                         0xEE  This type model the empty symbol
       COMPONENT_NULL                          0xFF Something was expected but nothing was found
       COMPONENT_ALGORITHM                     0xA0
       COMPONENT_INITIAL_SOLUTION_GENERATOR    0xB1
       COMPONENT_TERMINATION_CRITERION         0xB2
       COMPONENT_NEIGHBORHOOD                  0xB3
       COMPONENT_TABU_TENURE                   0xB4
       COMPONENT_NEIGHBORHOOD_OR_EMPTY         0xBE
       COMPONENT_PERTURBATION                  0xC1
       COMPONENT_ACCEPTANCE                    0xC2
       COMPONENT_SHAKE                         0xC3
       COMPONENT_NEIGHBORHOOD_CHANGE           0xC4
     */
    int type;
    /** @brief rawComponent
     * The pointer to the actual component
     */
    void* rawComponent;
    /** @brief token
     * The string token that represents this component (not used)
    */
    char* token;
public:
    /**
     * @brief Component
     *          Component constructor
     * @param type
     *          The type of the Component
     * @param rawData
     *          A void pointer that leads to an object of the type type
     */
    Component(int type,void* rawData):type(type),rawComponent(rawData) { }
    /**
     * @brief Component
     * The default component has type COMPONENT_NULL
     */
    Component():type(COMPONENT_NULL),rawComponent(nullptr) { }
    /**
     * @brief operator =
     *        Component copy operator
     * @param a
     *        The Component to copy
     * @return
     *        The copy.
     */
    virtual Component& operator=(const Component& a);    
    /**
     * @brief is
     *      Check if the component is of type "type".
     * @param type
     * @return
     * Return true if the component is of type "type".
     */
    bool is(int type) {return this->type == type;}
    /**
     * @brief getType
     *  returns Component type
     * @return
     */
    int getType() {return type;}
    /**
     * @brief setType
     *  sets Component type
     * @param type
     */
    void setType(int type) {this->type = type;}
    /**
     * @brief setRawComponent
     *  set the pointer to rawComponent
     * @param rawc
     */
    void setRawComponent(void* rawc) {this->rawComponent = rawc;}
    template<typename T>
    /**
     * @brief get
     * Casts rawCompoent to type T
     * @return
     */
    T* get(){ return (T*) rawComponent;}
};

class GeneralParserE;
/**
 * @brief The Builder class
 * A Builder is a class that is capable of building components for
 * one or more problems.
 * The Builder gets enrolled in the parsing process if it is compatible with a date problem
 *  and the GeneralParser will ask the builder to load the instance if it's able to do it.
 * All is required when extending this calss is the implementation of the method isCompatibleWith
 * and the overloading of all the buildX methods for which it's able to build Components.
 */
class Builder
{
protected:
    /**
     * @brief gp
     * The GeneralParserE instance
     */
    GeneralParserE& gp;
    /**
     * @brief tm
     * The TokenManager provided by GeneralParserE
     */
    TokenManager& tm;
    /**
     * @brief retrieveComponent
     *        The method to call to retrieve a component when needed.
     * @param type
     *         The type of the component to retrieve.
     *         See component documentation for the base types.
     * @param allowEmpty
     *          If set to true it could return an empty component otherwise it will
     *          rise an error.
     * @return
     *          A component of the type required ( if allowEmpty is true it could also
     *          be of type empty).
     */
    virtual Component retrieveComponent(int type,bool allowEmpty=false);
    /**
     * @brief buildNeighborhoodVector
     *          This utility method builds a vector of 1 to n Neighborhood
     *          reading from the token manager.
     * @return
     *          A vector of pointers to Neighborhood objects.
     */
    virtual std::vector<emili::Neighborhood*> buildNeighborhoodVector();
    /**
     * @brief buildComponentVector
     *          This utility method builds a vector of 1 to n Components
     *          reading from the token manager.
     * @return
     *          A vector of pointers to Components objects.
     */
    template <class T>
    std::vector<T*> buildComponentVector(int type);
public:
    /**
     * @brief Builder
     * @param generalParser
     *          The general Parser object that will use this builder.
     * @param tokenManager
     *          The tokenManager
     */
    Builder(GeneralParserE& generalParser,TokenManager& tokenManager):gp(generalParser),tm(tokenManager) { }
    /**
     * @brief openInstance
     *        If the Builder returns true to conOpenInstance, General parser can call this
     *        method to load the instance from a file.
     *        If this method is implemented in a new Builder,
     *        canOpenInstance(char *problem_definition) also has to be redefined otherwise
     *        GeneralParserE will never call the method.
     * @return
     *        A pointer to the problem instance object.
     */
    virtual emili::Problem* openInstance() {return nullptr;}
    /**
     * @brief isCompatibleWith
     *         This method tells GeneralParserE if this builder is compatible with a problem.
     *
     * @param problem_definition
     *          A string or char pointer that represents the problem definition.
     * @return
     *          true if the builder is compatible with the problem, false otherwise.
     */
    virtual bool isCompatibleWith(char* problem_definition)=0;
    /**
     * @brief isCompatibleWith
     *         This method tells GeneralParserE if this builder is compatible with a problem.
     *         The implementation of this method is mandatory if your are exenting this class.
     * @param problem_definition
     *          A string or char pointer that represents the problem definition.
     * @return
     *          true if the builder is compatible with the problem, false otherwise.
     */
    virtual bool isCompatibleWith(std::string& problem_definition) {return isCompatibleWith((char*)problem_definition.c_str());}
    /**
     * @brief canOpenInstance
     *         Tells GeneralParserE if the builder con load an instance of problem_definition
     * @param problem_definition
     *         A string or char pointer that represents the problem definition.
     * @return
     *         true if the builder can load the problem, false otherwise.
     */
    virtual bool canOpenInstance(char* problem_definition) {return false;}
    /**
     * @brief canOpenInstance
     *         Tells GeneralParserE if the builder con load an instance of problem_definition.
     *         This method has to be overridden when extending this class with the code to properly
     *         check when a new defined Builder is compatible with a problem.
     *
     * @param problem_definition
     *         A string or char pointer that represents the problem definition.
     * @return
     *         true if the builder can load the problem, false otherwise.
     */
    virtual bool canOpenInstance(std::string& problem_definition) {return canOpenInstance((char*)problem_definition.c_str());}
    /**
     * @brief typeName
     *         Returns a string that describe type.
     *         This method should be modified in
     *         a class that extends Builder and uses components different from the base ones
     * @param type
     * @return
     *        A string representing type.
     */
    virtual std::string typeName(int type){return std::string();}
    /**
     * @brief buildComponent
     *          Tries to build a component of type "type" using the token manager to
     *          retrieve the information about the specific class and parameters of the
     *          component.
     *          This method should be modified in a class that extends Builder and defines
     *          new component types.
     * @param type
     * @return
     *          A component of type "type" otherwise it returns a component of type COMPONENT_NULL
     *          if nothing is found or a component of type COMPONENT_EMPTY if the end of a rule
     *          that allows empty is reached.
     */
    virtual Component buildComponent(int type);
    /**
     * @brief buildAlgo
     *          This method is called by buildComponent(type) if type is COMPONENT_ALGORITHM.
     * @return
     *       a pointer to an object of type COMPONENT_ALGORITHM or, if nothing found, nullptr.
     */
    virtual emili::LocalSearch* buildAlgo() {return nullptr;}
    /**
     * @brief buildInitialSolution
     *          This method is called by buildComponent(type) if type is COMPONENT_INITIAL_SOLUTION_GENERATOR.
     * @return
     *       a pointer to an object of type COMPONENT_INITIAL_SOLUTION_GENERATOR or, if nothing found, nullptr.
     */
    virtual emili::InitialSolution* buildInitialSolution(){return nullptr;}
    /**
     * @brief buildAlgo
     *          This method is called by buildComponent(type) if type is COMPONENT_NEIGHBORHOOD.
     * @return
     *       a pointer to an object of type COMPONENT_NEIGHBORHOOD or, if nothing found, nullptr.
     *       At the end of a series of neighborhoods ( e.g. the specification of the neighborhoods
     *      used by a VNS) COMPONENT_EMPTY should be returned.
     */
    virtual emili::Neighborhood* buildNeighborhood(){return nullptr;}
    /**
     * @brief buildAlgo
     *          This method is called by buildComponent(type) if type is COMPONENT_TERMINATION_CRITERION.
     * @return
     *       a pointer to an object of type COMPONENT_TERMINATION_CRITERION or, if nothing found, nullptr.
     */
    virtual emili::Termination* buildTermination(){return nullptr;}
    /**
     * @brief buildAlgo
     *          This method is called by buildComponent(type) if type is COMPONENT_PERTURBATION.
     * @return
     *       a pointer to an object of type COMPONENT_PERTURBATION or, if nothing found, nullptr.
     */
    virtual emili::Perturbation* buildPerturbation(){return nullptr;}
    /**
     * @brief buildAcceptance
     *          This method is called by buildComponent(type) if type is COMPONENT_ACCEPTANCE.
     * @return
     *       a pointer to an object of type COMPONENT_ACCEPTANCE or, if nothing found, nullptr.
     */
    virtual emili::Acceptance* buildAcceptance(){return nullptr;}
    /**
     * @brief buildShake
     *          This method is called by buildComponent(type) if type is COMPONENT_SHAKE.
     * @return
     *       a pointer to an object of type COMPONENT_SHAKE or, if nothing found, nullptr.
     */
    virtual emili::Shake* buildShake(){return nullptr;}
    /**
     * @brief buildNeighborhoodChange()
     *          This method is called by buildComponent(type) if type is COMPONENT_NEIGHBORHOOD_CHANGE.
     * @return
     *       a pointer to an object of type COMPONENT_NEIGHBORHOOD_CHANGE or, if nothing found, nullptr.
     */
    virtual emili::NeighborhoodChange* buildNeighborhoodChange(){return nullptr;}
    /**
     * @brief buildAlgo
     *          This method is called by buildComponent(type) if type is COMPONENT_TABU_TENURE.
     * @return
     *       a pointer to an object of type COMPONENT_TABU_TENURE or, if nothing found, nullptr.
     */
    virtual emili::TabuMemory* buildTabuTenure(){return nullptr;}
    /**
     * @brief buildProblem
     *          This method is called by buildComponent(type) if type is COMPONENT_PROBLEM.
     * @return
     *       a pointer to an object of type COMPONENT_PROBLEM or, if nothing found, nullptr.
     */
    virtual emili::Problem* buildProblem(){return nullptr;}
    /**
     * @brief cmd_info
     *        This method should return a string representing all the supported components
     *        as well as the syntax to call them (e.g. additional literal parameters)
     *
     * @return
     *        info string.
     */
    virtual std::string cmd_info(){return std::string("");}
};
/**
 * @brief The GeneralParserE class
 *          This class defines the object responsible for the parsing of the command line arguments.
 *          It reads the instance file path, the problem definition and selects the compatible builders
 *          and guides the parsing process by calling the builders when needed. At the end it reads and
 *          set up the running time limit, the random seed and the print flag.
 */
class GeneralParserE: public GeneralParser
{
protected:
    /**
     * @brief allbuilders
     *         All the builder available at the beginning of the execution
     */
    std::vector< Builder* > allbuilders;
    /**
     * @brief activeBuilders
     *         This vector contains only the builders that are compatible with the problem definition.
     */
    std::vector< Builder* > activeBuilders;
    /**
     * @brief instance
     *         The problem instance
     */
    emili::Problem* instance;
    /**
     * @brief typeName
     *          Returns a string that describe type.
     * @param type
     * @return
     */
    virtual std::string typeName(int type);
public:
    /**
     * @brief GeneralParserE
     *         This constructor uses the parameters to initialize a TokenManager
     * @param tokens
     * @param numberOfTokens
     */
    GeneralParserE(char** tokens,int numberOfTokens):GeneralParser(tokens,numberOfTokens) { }
    GeneralParserE(TokenManager tokenmanager):GeneralParser(tokenmanager) { }
    /**
     * @brief parseParams
     *          This method parse, builds and returns the algorithm represented by the
     *          command line arguments used to call emili.
     * @return
     *          The algorithm ( if there is an error in the parsing the execution is terminated
     *           before the return)
     */
    virtual emili::LocalSearch* parseParams();
    /**
     * @brief addBuilder
     *          Adds b to the vector of available builders
     * @param b
     */
    virtual void addBuilder(Builder* b);
    /**
     * @brief getInstance
     *          returns a pointer to the problem instance object.
     * @return
     */
    virtual emili::Problem* getInstance() {return instance;}
    /**
     * @brief setInstance
     *          Sets instance to ins.
     * @param ins
     */
    virtual void setInstance(emili::Problem* ins) {instance=ins;}
    /**
     * @brief getTokenManager
     * @return
     *          returns a reference to the TokenManager.
     */
    virtual TokenManager& getTokenManager() {return tm;}
    /**
     * @brief buildComponent
     *         This method cycles trough activeBuilders until it finds
     *          one Builder that loads a compoent of type "type" using the
     *          TokenManager.
     * @param type
     * @return
     *         The component requested otherwise it generates an error and halts the execution.
     *
     */
    virtual Component buildComponent(int type);
    /**
     * @brief fatalError
     *          Prints of stdout an error occured during the parsing.
     * @param received_type
     * @param expected_type
     */
    virtual void fatalError(int received_type,int expected_type);    
};
/**
 * @brief The EmBaseBuilder class
 *     This class implements the Builder for the fundamental,
 *     problem indipendent components defined in emilibase.h.
 */
class EmBaseBuilder: public Builder
{
public:
    /**
     * @brief EmBaseBuilder
     * EmBaseBuilder loads the base components described in emilibase.h
     */
    EmBaseBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    /**
     * @brief isCompatibleWith
     * This Builder is compatible with any problem,
     * since all the components loaded are problem indipendent
     * @return
     *   true
     */
    virtual bool isCompatibleWith(char *problem_definition) { return true;}
    virtual emili::InitialSolution* buildInitialSolution();
    /**
     * @brief buildAlgo
     * This method will load emili::IteratedLocalSearch, emili::LocalSearch,
     * emili::VNDSearch, emili::BestTabuSearch, emili::FirstTabuSearch and trigger the load of all the other
     * components.
     * @return
     *  an algorithm of type LocalSearch
     */
    virtual emili::LocalSearch* buildAlgo();
    /**
     * @brief buildTermination
     * This method will load emili::LocalMinimaTermination, emili::MaxStepsTermination,
     * emili::WhileTrueTermination, emili::MaxStepsOrLocmin, emili::TimedTermination.
     * @return
     *  a termination
     */
    virtual emili::Termination* buildTermination();
    /**
     * @brief buildPerturbation
     * This method will load emili::RandomMovePerturbation, emili::VNRandomMovePerturbation and
     * emili::NoPerturbation.
     * @return
     *  a Perturbation
     */
    virtual emili::Perturbation* buildPerturbation();
    /**
     * @brief buildAcceptance
     * This method will load the acceptance criteria: emili::ImproveAccept, emili::AlwaysAccept,
     * emili::Metropolis, emili::AcceptPlateau.
     * @return
     */
    virtual emili::Acceptance* buildAcceptance();
    /**
     * @brief buildNeighborhood
     * This method will load the general neighborhood emili::RandomConstructiveHeuristicNeighborhood     *
     * @return
     * a neighborhood
     */
    virtual emili::Neighborhood* buildNeighborhood();
    /**
     * @brief buildShake
     * @return
     * a Shake operator to be used in a VNS algorithm
     */
    virtual emili::Shake* buildShake();
    /**
     * @brief buildNeighborhoodChange
     * @return
     * a NeighborhoodChange operator to be used in a VNS algorithm
     */
    virtual emili::NeighborhoodChange* buildNeighborhoodChange();
};

}

#endif
