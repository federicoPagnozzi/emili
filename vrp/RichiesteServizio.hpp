//
//  RichiesteServizio.hpp
//  OrCode
//
//  Created by garazi on 25/01/2017.
//  Copyright © 2017 Garazi. All rights reserved.
//

#ifndef RichiesteServizio_hpp
#define RichiesteServizio_hpp

#include <stdio.h>


class RichiesteServizio
{
public:
    int id; //identificativo richiesta di servizio
    int tipo; //1 PD-DO; 2 PO-DD; 3 PO-DO
    //string trattamento; //nome del trattamento
    int oc;// occupazione del veicolo (quanti posti vengono occupati)
    
    int timewinPmin; // minimo tempo Pick up
    int timewinPmax; // massimo tempo Pick up
    int timewinDmin; // minimo tempo Delivery
    int timewinDmax; // massimo tempo Delivery
    int locationP; //identificativo della localitá di pick up
    int locationD; //identificativo della localitá di delivery
    int paziente; //identificativo del paziente trasportato
    int idSosta; //se la sosta À di tipo A o di tipo B (1 o 2)
    int Barellato_SediaRotelle;//0 se non barellato (o sedia a rotelle) e 1 se barellato (o sedia a rotelle)
    int obbligo_associazione;//c'À il vincolo di dover soddisfare quella richieste solo con mezzi di una certa associazione? (1=si), (0=no)
    
    int staff;// occupazione del veicolo (quanti posti vengono occupati)
    int seated;// occupazione del veicolo (quanti posti vengono occupati)
    int stretcher;// occupazione del veicolo (quanti posti vengono occupati)
    int wheelchair;// occupazione del veicolo (quanti posti vengono occupati)
    int st;
    double RideTime;
    
    
public:
    
    void display_ric();
};






#endif /* RichiesteServizio_hpp */
