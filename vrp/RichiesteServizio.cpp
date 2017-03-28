//
//  RichiesteServizio.cpp
//  C-B
//
//  Created by garazi on 01/02/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include "RichiesteServizio.hpp"

#include<iostream>

void RichiesteServizio::display_ric(){
    
    std::cout << id << " " << tipo << " " << timewinPmin << " " << timewinPmax << " " << timewinDmin << " " << timewinDmax << " " << locationP << " " << locationD << " " << staff << " " << seated << " " << stretcher << " " << wheelchair << " " << st << " " << RideTime << std::endl;
    
}

