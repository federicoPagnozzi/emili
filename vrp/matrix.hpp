//
//  matrix.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <stdio.h>

#include"RichiesteServizio.hpp"
#include <vector>
int* timewinminP_order(std::vector<RichiesteServizio*> ric, int n, int*E);
int* random_order(std::vector<RichiesteServizio*> ric, int numRichieste, int*E);
int* random_order2(int m, int*E);
#endif /* matrix_hpp */