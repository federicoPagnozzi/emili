//
//  Veicoli.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef Veicoli_hpp
#define Veicoli_hpp

#include <stdio.h>


class Veicoli
{
public:
    int id; //identificativo del veicolo
    int location; //identificativo della localit‡ del deposito del veicolo
    int cap; // capacit‡ del veicolo
    int disponibili; //numero di veicoli disponibili
    int tipo; // tipologia a cui appartiene il veicolo
    double distanza; //distanza oltre la quale non si ha pi˘ un costo fisso, ma un costo variabile in base ai km percorsi
    double costo1; //Costo del veicolo fisso fino a "distanza" chilometri
    double costo2; //Costo variabile (al km) dopo "distanza" chilometri
    double costoNumeroPazientiTrasportati;
    int staff; // capacit‡ del veicolo
    int seated; // capacit‡ del veicolo
    int stretcher; // capacit‡ del veicolo
    int wheelchair; // capacit‡ del veicolo
    double MaxRoute;
    
    
public:
    
    void display_veic();
    
}; // fine classe veicoli




#endif /* Veicoli_hpp */
