//
//  openfiles.cpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#include "openfiles.hpp"

#include<iostream>
#include<fstream>
#include<string>
#include <sstream>
#include<cmath>
#include <cstring>
#include <ctime>
//#include <sys/time.h>
#include "openfiles.hpp"
#include "RichiesteServizio.hpp"
#include "Veicoli.hpp"

using namespace std;

double** openf(int m,int n, string a, string b,double**A)
{
    //double** A=0;
    //A=new double*[m];
    
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+b+"/"+a+".txt";
    ifstream ar1(c.c_str());
	   int i,j;
    
	   while (!ar1.eof()) {
           for(i=0; i<m; i++)
           { //A[i]=new double[n];
               for(j=0; j<n; j++){
                   ar1 >> A[i][j];}
           }} //A matrix is the matrix of dimensions
    /*
     for(i=0; i<m; i++){
     for(j=0; j<n; j++){cout << A[i][j] << " " ;
     }cout << endl;}*/
	   ar1.close();
    return A;
}
int** openi(int m,int n, string a, string b,int**AA)
{
    //int** AA=0;
    //AA=new int*[m];
    
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+b+"/"+a+".txt";
    ifstream ar11(c.c_str());
	   int i,j;
    
	   while (!ar11.eof()) {
           for(i=0; i<m; i++)
           { //AA[i]=new int[n];
               for(j=0; j<n; j++){
                   ar11 >> AA[i][j];}
           }} //A matrix is the matrix of dimensions
    /*for(i=0; i<m; i++){
     for(j=0; j<n; j++){cout << AA[i][j] << " " ;
     }cout << endl;}*/
	   ar11.close();
    return AA;
}

double opencost(int m, int n, string a, string b,double cost)
{
    double** A;
    //AA=new int*[m];
    A=new double*[m];
    int i,j;
    for (i=0; i<m; i++)
    {A[i]=new double[1];}
    
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+b+"/"+a+".txt";
    ifstream ar11(c.c_str());
    
	   for(i=0; i<m; i++)
       {A[i][0]=0;}
    
	   while (!ar11.eof()) {
           for(i=0; i<m; i++)
           { for(j=0; j<1;j++){
               ar11 >> A[i][j];}}
       } //A matrix is the matrix of dimensions
    /*for(i=0; i<m; i++){
     for(j=0; j<n; j++){cout << AA[i][j] << " " ;
     }cout << endl;}*/
    
	   cost=A[n][0];
    
	   ar11.close();
    
	   for(i=0; i<m; i++)
       {delete[] A[i];}
	   delete[] A;
    
    return cost;
}
int* openv(int m, string a, string b,int*AA)
{
    //int** AA=0;
    //AA=new int*[m];
    
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+b+"/"+a+".txt";
    ifstream ar11(c.c_str());
	   int i;
    
	   while (!ar11.eof()) {
           for(i=0; i<m; i++)
           {
               ar11 >> AA[i];
           }}
	   ar11.close();
    return AA;
}
double*** opens(int m,int n, int l, string a, string b,double***A)
{
    //double** A=0;
    //A=new double*[m];
    
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+b+"/"+a+".txt";
    ifstream ar1(c.c_str());
	   int i,j,k;
    
	   while (!ar1.eof()) {
           for(i=0; i<m; i++)
           { //A[i]=new double[n];
               for(k=0; k<l; k++){
                   for(j=0; j<n; j++){
                       ar1 >> A[i][j][k];}
               }}} //A matrix is the matrix of dimensions
    /*
     for(i=0; i<m; i++){
     for(j=0; j<n; j++){cout << A[i][j] << " " ;
     }cout << endl;}*/
	   ar1.close();
    return A;
}


void createsolfile(double totalcost, int totalvehicles, int numVeicoli, int numRichieste, double** Solreq, double** ind, double** E, double** C, string a)
{
    
    
    ofstream ar15;
    string b=a+".txt";
    ar15.open(b.c_str());
    int i, j;
    ar15 << "COSTE:  " << totalcost << endl;
    ar15 << "VEHICLES USED:  " << totalvehicles << endl << endl;
    for (i=0; i<numVeicoli; i++)
    { ar15 << "Ruta relativa al vehiculo   " << i << endl;
        ar15 << endl;
        for (j=1; j<2*numRichieste+2; j++)
        {
            if(Solreq[i][j]==-2){}else{
                
                if (Solreq[i][j]==-1){
                    ar15 << 2*numRichieste+i << " " << 0 << " " << 0 << " " << 0 <<" " << 24*60 << " " << C[i][1] << " " << -1 << " "<< -1 << " " << -1;
                    
                }else{ if(ind[i][j]==1){
                    ar15 << Solreq[i][j] << " ";
                    if(E[int(Solreq[i][j])][1]==1){ar15 << 1 << " ";}
                    else{ar15 << 3 << " ";}
                    ar15 <<  1 << " " << E[int(Solreq[i][j])][3] << " " << E[int(Solreq[i][j])][4] << " " << E[int(Solreq[i][j])][7] << " " << E[int(Solreq[i][j])][9] << " " << E[int(Solreq[i][j])][10] << " " << E[int(Solreq[i][j])][11];
                }else{
                    ar15 << numRichieste+Solreq[i][j] << " ";
                    if(E[int(Solreq[i][j])][1]==1){ar15 << 2 << " ";}else{ar15 << 4 << " ";}
                    ar15 << -1 << " " << E[int(Solreq[i][j])][5] << " " << E[int(Solreq[i][j])][6] << " " << E[int(Solreq[i][j])][8] << " " << E[int(Solreq[i][j])][9] << " " << E[int(Solreq[i][j])][10] << " " << E[int(Solreq[i][j])][11];
                }}} ar15 << "\t\t" ;
        }
        ar15 << endl;
    }
    ar15.close();
    
}

void guardarmatriz(double*** Sol, int numVeicoli, int numRichieste, double cost, double f_cost, int sum, string w, int contt, string name, string h)
{
    ofstream ar16;
    
    
    
    
    string b="C:/Users/pc/Desktop/VNS/DATA_newCost/"+h+"/"+w+"solucion nueva"+name+".txt";
    
    //b << "solucion nueva.txt" ;
    
    
    
    ar16.open(b.c_str() );
    
    
    ar16 << "REAL COST WITH VIOLATION-->  " << cost << endl << endl << endl;
    ar16 << "COST WITH VIOLATIONS -->" << f_cost << endl << endl;
    ar16 << "NUMBER OF VEHICLES USED " << sum << endl;
    ar16 << "This solution was found on the iteration no. " << contt << endl;
    
    for (int i=0; i<numVeicoli; i++){
        
        ar16 << "RUTA " << i << endl;
        for(int j=0; j<2*numRichieste; j++)
        {for (int z=0; z<4; z++){
            if(Sol[i][j][z]!=-2){
                ar16 << Sol[i][j][z] << "  " ;
            }}if(Sol[i][j][0]!=-2){ar16 << endl;
            }
        }
        
        
    }
    ar16.close();
    return; }

void guardarmatrizadjusted(double*** Sol, int numVeicoli, int numRichieste, double cost, double f_cost, int sum, string w, int contt, string name, double f_costbest, string h)
{
    ofstream ar16;
    
    
    
    string b="C:/Users/pc/Desktop/VNS/DATA_newCost/"+h+"/"+w+"solucion nueva"+name+".txt";
    
    //b << "solucion nueva.txt" ;
    
    ar16.open(b.c_str() );
    
    ar16<< "This is the adjusted solution taken from the solution on iteration no. " << contt << endl;
    ar16 << "The previous solution had a cost of " << f_costbest << endl;
    ar16 << "REAL COST WITH VIOLATION-->  " << cost << endl << endl << endl;
    ar16 << "COST WITH VIOLATIONS -->" << f_cost << endl << endl;
    ar16 << "NUMBER OF VEHICLES USED " << sum << endl;
    
    
    for (int i=0; i<numVeicoli; i++){
        
        ar16 << "RUTA " << i << endl;
        for(int j=0; j<2*numRichieste; j++)
        {for (int z=0; z<4; z++){
            if(Sol[i][j][z]!=-2){
                ar16 << Sol[i][j][z] << "  " ;
            }}if(Sol[i][j][0]!=-2){ar16 << endl;
            }
        }
        
        
    }
    ar16.close();
    return; }

void guardarmatrizI(double*** Sol, int numVeicoli, int numRichieste, double cost, int sum)
{
    ofstream ar16;
    
    
    
    stringstream b;
    
    b << "solucion inicial.txt" ;
    string g=b.str();
    ar16.open(g.c_str() );
    
    
    ar16 << "COSTE -->  " << cost << endl << endl << endl;
    ar16 << "NUMBER OF VEHICLES USED -->" << sum << endl << endl;
    for (int i=0; i<numVeicoli; i++){
        
        ar16 << "RUTA " << i << endl;
        for(int j=0; j<2*numRichieste; j++)
        {for (int z=0; z<4; z++){
            if(Sol[i][j][z]!=-2){
                ar16 << Sol[i][j][z] << "  " ;
            }}if(Sol[i][j][0]!=-2){ar16 << endl;
            }
        }
        
        
    }
    ar16.close();
    return; }

RichiesteServizio* openric(int m,string a, string b,RichiesteServizio* ric){
    
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+b+"/"+a+".txt";
    ifstream ar11(c.c_str());
    int i;
    
    
    
    while (!ar11.eof()) {
        for(i=0; i<m; i++)
        {
            
            ar11 >> ric[i].id;
            ar11 >> ric[i].tipo;
            ar11 >> ric[i].oc;
            ar11 >> ric[i].timewinPmin;
            ar11 >> ric[i].timewinPmax;
            ar11 >> ric[i].timewinDmin;
            ar11 >> ric[i].timewinDmax;
            ar11 >> ric[i].locationP;
            ar11 >> ric[i].locationD;
            ar11 >> ric[i].paziente;
            ar11 >> ric[i].idSosta;
            ar11 >> ric[i].Barellato_SediaRotelle;
            ar11 >> ric[i].obbligo_associazione;
            ar11 >> ric[i].staff;
            ar11 >> ric[i].seated;
            ar11 >> ric[i].stretcher;
            ar11 >> ric[i].wheelchair;
        }} //A matrix is the matrix of dimensions
    /*
     for(i=0; i<m; i++){
     for(j=0; j<n; j++){cout << A[i][j] << " " ;
     }cout << endl;}*/
    ar11.close();
    
    return ric;
}

Veicoli* openveic(int m,string a,string b, Veicoli* veic){
    string c="/Users/garazi/Documents/Garazi/Data/new_instances/"+b+"/"+a+".txt";
    ifstream ar11(c.c_str());
    int i;
    
    
    
    while (!ar11.eof()) {
			   		for(i=0; i<m; i++)
                    {
                        
                        
                        ar11 >> veic[i].id;
                        ar11 >> veic[i].location;
                        ar11 >> veic[i].cap;
                        ar11 >> veic[i].disponibili;
                        ar11 >> veic[i].tipo;
                        ar11 >> veic[i].distanza;
                        ar11 >> veic[i].costo1;
                        ar11 >> veic[i].costo2;
                        ar11 >> veic[i].costoNumeroPazientiTrasportati;
                        ar11 >> veic[i].staff;
                        ar11 >> veic[i].seated;
                        ar11 >> veic[i].stretcher;
                        ar11 >> veic[i].wheelchair;
                        
                    }} //A matrix is the matrix of dimensions
    /*
     for(i=0; i<m; i++){
     for(j=0; j<n; j++){cout << A[i][j] << " " ;
     }cout << endl;}*/
    ar11.close();
    
    
    
    
    return veic;
}
