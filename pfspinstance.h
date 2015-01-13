/***************************************************************************
 *   Copyright (C) 2012 by Jérémie Dubois-Lacoste   *
 *   jeremie.dl@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef _PFSPINSTANCEWT_H_
#define _PFSPINSTANCEWT_H_

#include <string>
#include <vector>

using namespace std;

class PfspInstance{
  private:
    int nbJob;
    int nbMac;
    std::vector< long int > dueDates;
    std::vector< long int > priority;
    bool silence;
    std::vector< std::vector <long int> > processingTimesMatrix;

  public:
    PfspInstance();
    ~PfspInstance();

    /* Read write privates attributs : */
    int getNbJob();
    int getNbMac();

    /* Allow the memory for the processing times matrix : */
    void allowMatrixMemory(int nbJ, int nbM);

    /* Read\Write values in the matrix : */
    long int getTime(int job, int machine);
    void setTime(int job, int machine, long int processTime);

    long int getDueDate(int job);
    void setDueDate(int job, int value);

    long int getPriority(int job);
    void setPriority(int job, int value);

    /* Read Data from a file : */
    bool readDataFromFile(char * fileName);
    /*Compute weighted tardiness*/
    long int computeWT (vector< int > & sol);
    /*Compute partial weighted tardiness*/
    long int computeWT (vector< int > & sol, int size);
    /* Compute MakeSpan */
    long int computeMS (vector<int> & sol);

    /* Compute partial MakeSpan*/
    long int computeMS(vector<int> &sol,int size);

    void setSilence(bool s);

};

#endif
