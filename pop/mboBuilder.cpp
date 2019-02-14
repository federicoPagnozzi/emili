//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#include "mboBuilder.h"
#include "vig_de.h"
#include "maneh.h"
#include "hdtlm.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <iostream>
#include <sstream>
#include <algorithm>

/* Algos */
#define LS_EMBO "embo"
#define LS_VigDE "vigde"
#define LS_MANEH "maneh"

emili::LocalSearch* prs::MboBuilder::buildAlgo()
{
    prs::incrementTabLevel();
    emili::LocalSearch* ls = nullptr;
    if(tm.checkToken(LS_EMBO))
    {
        printTab("EMBO");

        emili::LocalSearch* ls1 = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::Neighborhood* ne = retrieveComponent(COMPONENT_NEIGHBORHOOD).get<emili::Neighborhood>();        
        emili::Perturbation* p = retrieveComponent(COMPONENT_PERTURBATION).get<emili::Perturbation>();
        int popsize = tm.getInteger();
        printTabPlusOne("Popsize",popsize);
        int kp = tm.getInteger();
        printTabPlusOne("k",kp);
        int xp = tm.getInteger();
        printTabPlusOne("x",xp);
        int mp = tm.getInteger();
        printTabPlusOne("m",mp);
        int age = tm.getInteger();
        printTabPlusOne("age",age);
        int q0 = tm.getDecimal();
        printTabPlusOne("q0",q0);
        ls = new emili::pop::EMBO(*in,*te,*ne,p,ls1,popsize,kp,xp,mp,age,q0);        
    }
    else if(tm.checkToken(LS_VigDE))
    {
        printTab("VigDE");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::InitialSolution* in1 = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        int popsize = tm.getInteger();
        printTabPlusOne("Popsize",popsize);
        float mf = tm.getDecimal();
        printTabPlusOne("Mutation factor",mf);
        float cr = tm.getDecimal();
        printTabPlusOne("Crossover factor",cr);
        ls = new emili::pop::vIG_DE(*te,*in,*in1,popsize,mf,cr);
    }
    else if(tm.checkToken(LS_MANEH))
    {
        printTab("MANEH");
        emili::InitialSolution* in = retrieveComponent(COMPONENT_INITIAL_SOLUTION_GENERATOR).get<emili::InitialSolution>();
        emili::Termination* te = retrieveComponent(COMPONENT_TERMINATION_CRITERION).get<emili::Termination>();
        emili::LocalSearch* ls1 = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        emili::LocalSearch* ls2 = retrieveComponent(COMPONENT_ALGORITHM).get<emili::LocalSearch>();
        emili::Perturbation* pe = retrieveComponent(COMPONENT_PERTURBATION).get<emili::Perturbation>();
        int pops = tm.getInteger();
        printTabPlusOne("Population size",pops);
        float pmp = tm.getDecimal();
        printTabPlusOne("Mutation probability",pmp);
        float pcp = tm.getDecimal();
        printTabPlusOne("Crossover probability",pcp);
        float bratio = tm.getDecimal();
        printTabPlusOne("bratio",bratio);
        float lsp = tm.getDecimal();
        printTabPlusOne("Local search probability",lsp);
        float alph = tm.getDecimal();
        printTabPlusOne("Alpha",alph);
        ls = new emili::pop::MANEH(*te,*in,pops,ls1,ls2,pe,pmp,pcp,bratio,lsp,alph);
    }


    prs::decrementTabLevel();
    return ls;
}

/*
extern "C" {
    prs::Builder* getBuilder(prs::GeneralParserE* ge)
    {
        return new prs::MboBuilder(*ge,ge->getTokenManager());
    }
}
*/
