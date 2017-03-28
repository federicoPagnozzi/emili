//
//  openfiles.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef openfiles_hpp
#define openfiles_hpp
#include <string>
#include <stdio.h>
#include "RichiesteServizio.hpp"
#include "Veicoli.hpp"
using namespace std;

double** openf(int m,int n,string a, string b,double**A);
int** openi(int m,int n, string a, string b,int**AA);
int* openv(int m, string a, string b,int*AA);
double opencost(int m, int n, string a, string b,double cost);
double*** opens(int m,int n, int l, string a, string b,double***A);
void createsolfile(double totalcost, int totalvehicles, int numVeicoli, int numRichieste, double** Solreq, double** ind, double** E, double** C, string a);
void guardarmatriz(double*** Sol, int numVeicoli, int numRichieste, double cost, double f_cost,int sum, string w, int contt, string name, string h);
void guardarmatrizI(double*** Sol, int numVeicoli, int numRichieste, double cost, int sum);
void guardarmatrizadjusted(double*** Sol, int numVeicoli, int numRichieste, double cost, double f_cost, int sum, string w, int contt, string name, double f_costbest, string h);
RichiesteServizio* openric(int m, string a, string b, RichiesteServizio* ric);
Veicoli* openveic(int m,string a,string b, Veicoli* veic);
#endif /* openfiles_hpp */
