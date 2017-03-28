//
//  actualizar.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include "actualizar.hpp"

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include <cstring>
#include <ctime>
#include <float.h>
#include <cfloat>

#include "openfiles.hpp"
#include "matrix.hpp"


#include "RichiesteServizio.hpp"




RichiesteServizio* TimeWindowTightening(int numRichieste, RichiesteServizio* ric, double MaxRideTime, double MinAt, double** MatTemp){
    
    int i;
    for(i=0; i<numRichieste;i++){
        if(ric[i].tipo==1){
            if(ric[i].timewinDmin-MaxRideTime-MinAt<0){
                ric[i].timewinPmin=0;
            }else{
                ric[i].timewinPmin=ric[i].timewinDmin-MaxRideTime-MinAt;
            }
            
            if(ric[i].timewinDmax-MinAt-MatTemp[ric[i].locationP][ric[i].locationD]>1440){
                ric[i].timewinPmax=1440;
            }else{
                ric[i].timewinPmax=ric[i].timewinDmax-MinAt-MatTemp[ric[i].locationP][ric[i].locationD];
            }
        }
        if(ric[i].tipo==2){
            if(ric[i].timewinPmin+MinAt+MatTemp[ric[i].locationP][ric[i].locationD]<0){
                ric[i].timewinDmin=0;
            }else{
                ric[i].timewinDmin=ric[i].timewinPmin+MinAt+MatTemp[ric[i].locationP][ric[i].locationD];
            }
            if(ric[i].timewinPmax+MinAt+MaxRideTime>1440){
                ric[i].timewinDmax=1440;
            }else{
                
                ric[i].timewinDmax=ric[i].timewinPmax+MinAt+MaxRideTime;
            }
            
            
            
        }
        
        
        
    }
    
    
    return ric;
}

