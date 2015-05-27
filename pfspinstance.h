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
    void computeTails(std::vector<int> &sol, int size,std::vector< std::vector< int > > & tail);
  public:
    PfspInstance(PfspInstance& is)
    {
        this->nbJob = is.getNbJob();
        this->nbMac = is.getNbMac();
        this->dueDates = is.getDueDates();
        this->priority = is.getPriorities();
        this->processingTimesMatrix = is.getProcessingTimesMatrix();
        this->silence = is.silence;
    }

    PfspInstance();
    ~PfspInstance();

    /* Read write privates attributes : */
    int getNbJob();
    void setNbJob(int jobCount);

    int getNbMac();
    void setNbMac(int machienCount);

    /* Allow the memory for the processing times matrix : */
    void allowMatrixMemory(int nbJ, int nbM);

    /* Read\Write values in the matrix : */
    long int getTime(int job, int machine);
    void setTime(int job, int machine, long int processTime);

    long int getDueDate(int job);
    void setDueDate(int job, int value);

    long int getPriority(int job);
    void setPriority(int job, int value);

    std::vector< long int >& getDueDates()
    {
        return dueDates;
    }

    std::vector< long int >& getPriorities()
    {
        return priority;
    }

    /* Read Data from a file : */
    bool readDataFromFile(char * fileName);

    /* Read Data from file with the other format*/
    bool readDataFromFile(const std::string _filename);

    /*Compute weighted tardiness*/
    long int computeWT (vector< int > & sol);
    /*Compute partial weighted tardiness*/
    long int computeWT (vector< int > & sol, int size);
    /* Compute MakeSpan */
    long int computeMS (vector<int> & sol);
    /* Compute FlowTime*/
    long int computeFT(vector<int> & sol);
    long int computeFT(vector<int> &sol, int size);
    /*Compute weighted completion times*/    
    long int computeWCT (vector< int > & sol);
    /*Compute partial weighted completion times*/    
    long int computeWCT (vector< int > & sol, int size);
    /* Compute partial MakeSpan*/
    long int computeMS(vector<int> &sol,int size);
    /*Compute weighted earliness*/
    long int computeWE (vector< int > & sol);
    /*Compute partial weighted earliness*/
    long int computeWE (vector< int > & sol, int size);
    /*Compute tardiness*/
    long int computeT(vector< int > & sol);
    /*Compute partial tardiness*/
    long int computeT(vector< int > & sol, int size);
    /*Compute earliness*/
    long int computeE (vector< int > & sol);
    /*Compute partial earliness*/
    long int computeE (vector< int > & sol, int size);
    /*Compute total completion time*/
    long int computeTCT(vector< int > &sol);
    long int computeTCT(vector< int > &sol,int size);

    /*Compute no wait make span*/
    long int computeNWMS(vector< int > & sol);
    /*Compute partial no wait make span*/
    long int computeNWMS(vector<int> & sol, int size);
    /*Compute no wait weighted tardiness*/
    long int computeNWWT(vector<int> &sol);
    long int computeNWWT(vector<int> &sol,int size);
    /*Compute no wait weighted earliness*/
    long int computeNWWE(vector<int> &sol);
    long int computeNWWE(vector<int> &sol,int size);
    /*Compute no wait earliness*/
    long int computeNWE(vector<int> &sol);
    long int computeNWE(vector<int> &sol,int size);
    /*Compute no wait tardiness*/
    long int computeNWT(vector<int> &sol);
    long int computeNWT(vector<int> &sol,int size);
    /*Compute no wait total completion time*/
    long int computeNWTCT(vector< int > &sol);
    long int computeNWTCT(vector< int > &sol,int size);

    /*Compute no idle make span*/
    long int computeNIMS(vector<int> & sol);
    /*Compute no idle partial make span*/
    long int computeNIMS(vector<int> &sol, int size);
    /*Compute no idle make span iwthout computing the sums of machine 1 processing times*/
    long int computeNIMS(vector<int> &sol, long int nims);
    /*Compute no idle weighted tardiness*/
    long int computeNIWT(vector<int> &sol);
    long int computeNIWT(vector<int> &sol,int size);
    /*Compute no idle weighted earliness*/
    long int computeNIWE(vector<int> &sol);
    long int computeNIWE(vector<int> &sol,int size);
    /*Compute no idle earliness*/
    long int computeNIE(vector<int> &sol);
    long int computeNIE(vector<int> &sol,int size);
    /*Compute no idle tardiness*/
    long int computeNIT(vector<int> &sol);
    long int computeNIT(vector<int> &sol,int size);
    /*Compute no idle total completion time*/
    long int computeNITCT(vector< int > &sol);
    long int computeNITCT(vector< int > &sol,int size);

    /* Compute weighted tardines starting from an index*/
    long int computeWT(vector< int > & sol, vector< vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i);

    long int computeWT(vector<int> &sol,vector<int>& prevJob,int job,vector<int>& previousMachineEndTime);
    //to compute WT and save the values for the job job
    void computeWTs(vector<int> &sol,vector<int>& prevJob,int job,vector<int>& previousMachineEndTime);
    void setSilence(bool s);

    //Compute the earliest completion time starting from the beginning and starting from the end in order to apply taillard acceleration...
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    void computeTails(std::vector<int> &sol, std::vector < std::vector< std::vector< int > > >& tails);
    void computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    int computeIdleTimeCoeff(vector<int>& prevJob, int job);

    const std::vector< std::vector < long int > > & getProcessingTimesMatrix() { return processingTimesMatrix; }



};

#endif
