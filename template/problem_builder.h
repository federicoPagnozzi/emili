 //
//  Created by Federico Pagnozzi on 12/12/17.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef ProblemXBUILDER_H
#define  ProblemXBUILDER_H
#include "../generalParser.h"
namespace prs
{
namespace problemX
{
/**
 * @brief The ProblemXBuilder class
 * A Builder is a class that is capable of building components for
 * one or more problems.
 * The Builder gets enrolled in the parsing process if it is compatible with a date problem
 *  and the GeneralParser will ask the builder to load the instance if it's able to do it.
 * All is required when extending this calss is the implementation of the method isCompatibleWith
 * and the overloading of all the buildX methods for which it's able to build Components.
 * This class models a builder for ProblemX and all the components implemented for this problem
 * In this template we override only the methods to create an InitialSolution, a Neighborhood and a Perturbation but,
 * methods for other components can be created as needed following the methods in the Builder class.
 * To activate the builder, a builder object must be created and added to the builders pool of GeneralParser in main.cpp
 * For an example of how to do it properly see the commented code at line 98 and 100 of main.cpp where ProblmXBuilder is activated
 */
class ProblemXBuilder: public Builder
{
public:
    /**
     * @brief ProblemXBuilderBuilder
     * @param generalParser
     *          The general Parser object that will use this builder.
     * @param tokenManager
     *          The tokenManager
     */
    ProblemXBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):Builder(generalParser,tokenManager) { }
    /**
     * @brief isCompatibleWith
     *         This method tells GeneralParserE if this builder is compatible with a problem.
     *
     * @param problem_definition
     *          A string or char pointer that represents the problem definition.
     * @return
     *          true if the builder is compatible with the problem, false otherwise.
     */
    virtual bool isCompatibleWith(char* problem_definition);
    /**
     * @brief canOpenInstance
     *         Tells GeneralParserE if the builder con load an instance of problem_definition
     * @param problem_definition
     *         A string or char pointer that represents the problem definition.
     * @return
     *         This method shoudl return true if the problem_definition corresponds to problemX
     *         and false otherwise.
     */
    virtual bool canOpenInstance(char* problem_definition);
    /**
     * @brief openInstance
     *       General parser can call this
     *        method to load the instance of ProblemX specified in the parameters.
     *        If this method is implemented in a new Builder,
     *        canOpenInstance(char *problem_definition) also has to be redefined otherwise
     *        GeneralParserE will never call the method.
     * @return
     *        A pointer to the problem instance object.
     */
    virtual emili::Problem* openInstance();

    virtual emili::InitialSolution* buildInitialSolution();

    virtual emili::Neighborhood* buildNeighborhood();

    virtual emili::Perturbation* buildPerturbation();
};
}
}
#endif //  PFSPBUILDER_H
