//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

#ifndef MBOBUILDER_H
#define  MBOBUILDER_H
#include "../generalParser.h"
#include "mbo.h"
namespace prs
{

class MboBuilder: public Builder
{
public:
    MboBuilder(GeneralParserE& generalParser,TokenManager& tokenManager):
        Builder(generalParser,tokenManager) { }
    virtual bool isCompatibleWith(char* problem_definition) {return true;}
    virtual bool canOpenInstance(char* problem_definition) {return false;}
    virtual emili::LocalSearch* buildAlgo();    
};
}
#endif //  MBOBUILDER_H
