//
//  actualizar.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef actualizar_hpp
#define actualizar_hpp

#include <stdio.h>
#include "RichiesteServizio.hpp"


RichiesteServizio* TimeWindowTightening(int numRichieste, RichiesteServizio* ric, double MaxRideTime, double MinAt, double** MatTemp);
#endif /* actualizar_hpp */
