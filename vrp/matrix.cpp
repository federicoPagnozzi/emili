//
//  matrix.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include "matrix.hpp"


#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include <cstring>
#include <ctime>
#include"RichiesteServizio.hpp"
#include"../emilibase.h"
//#include <sys/time.h>
//#include "openfiles.h"


using namespace std;

int* timewinminP_order(std::vector<RichiesteServizio*> ric, int n, int*E)
{
    int i, j;
    
    int aux=0;
    //int* E;
    int *G;
    //E=new int[n];
    G=new int[n];
    
    for(i=0; i<n; i++)
    {G[i]=ric[i]->timewinPmin;
        E[i]=ric[i]->id;
    }
    
    for(i=0; i<n-1; i++)
    {
        for (j=i+1; j<n; j++)
        {
            if(G[i]>G[j])
            {
                aux=E[i];
                E[i]=E[j];
                E[j]=aux;
                aux=G[i];
                G[i]=G[j];
                G[j]=aux;
            }
        }
    }
    
    delete [] G;
    
    return E;
}


int* random_order(std::vector<RichiesteServizio*> ric, int numRichieste, int*E)
{
    int i;
    
    int aux=0;
    
    
    
    for(i=0; i<numRichieste; i++)
    {E[i]=ric[i]->id;
    }
    
    for(i=0; i<numRichieste; i++)
    {
        //int r=rand() % numRichieste;
        int r=emili::generateRandomNumber() % numRichieste;
        aux=E[i];
        E[i]=E[r];
        E[r]=aux;
    }
    
    
    
    return E;
}
int* random_order2(int m, int*E){
    int i;
    
    int aux=0;
    
    
    
    for(i=0; i<m; i++)
    {
        int r=emili::generateRandomNumber() % m;
        aux=E[i];
        E[i]=E[r];
        E[r]=aux;
    }
    
    
    
    return E;
}
