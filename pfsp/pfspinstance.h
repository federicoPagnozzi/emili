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



class PfspInstance{
  private:
    int nbJob;
    int nbMac;
    std::vector< long int > dueDates;
    std::vector< long int > priority;
    bool silence;
    std::vector< std::vector <long int> > processingTimesMatrix;
    std::vector< std::vector < std::vector< int > > > setUpTimes;
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
        this->setUpTimes = is.setUpTimes;
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

    /* Read Data from sequence dependent setup times file*/
    bool readSeqDepDataFromFile(char* filename);

    /*Compute weighted tardiness*/
    long int computeWT (std::vector< int > & sol);
    /*Compute partial weighted tardiness*/
    long int computeWT (std::vector< int > & sol, int size);
    /* Compute MakeSpan */
    long int computeMS (std::vector<int> & sol);
    /* Compute FlowTime*/
    long int computeFT(std::vector<int> & sol);
    long int computeFT(std::vector<int> &sol, int size);
    /*Compute weighted completion times*/    
    long int computeWCT (std::vector< int > & sol);
    /*Compute partial weighted completion times*/    
    long int computeWCT (std::vector< int > & sol, int size);
    /* Compute partial MakeSpan*/
    long int computeMS(std::vector<int> &sol,int size);
    /*Compute weighted earliness*/
    long int computeWE (std::vector< int > & sol);
    /*Compute partial weighted earliness*/
    long int computeWE (std::vector< int > & sol, int size);
    /*Compute tardiness*/
    long int computeT(std::vector< int > & sol);
    /*Compute partial tardiness*/
    long int computeT(std::vector< int > & sol, int size);
    /*Compute earliness*/
    long int computeE (std::vector< int > & sol);
    /*Compute partial earliness*/
    long int computeE (std::vector< int > & sol, int size);
    /*Compute total completion time*/
    long int computeTCT(std::vector< int > &sol);
    long int computeTCT(std::vector< int > &sol,int size);

    /*Compute no wait make span*/
    long int computeNWMS(std::vector< int > & sol);
    /*Compute partial no wait make span*/
    long int computeNWMS(std::vector<int> & sol, int size);
    /*Compute no wait weighted tardiness*/
    long int computeNWWT(std::vector<int> &sol);
    long int computeNWWT(std::vector<int> &sol,int size);
    /*Compute no wait weighted earliness*/
    long int computeNWWE(std::vector<int> &sol);
    long int computeNWWE(std::vector<int> &sol,int size);
    /*Compute no wait earliness*/
    long int computeNWE(std::vector<int> &sol);
    long int computeNWE(std::vector<int> &sol,int size);
    /*Compute no wait tardiness*/
    long int computeNWT(std::vector<int> &sol);
    long int computeNWT(std::vector<int> &sol,int size);
    /*Compute no wait total completion time*/
    long int computeNWTCT(std::vector< int > &sol);
    long int computeNWTCT(std::vector< int > &sol,int size);
    long int computeNWWCT(std::vector< int > &sol);
    long int computeNWWCT(std::vector< int > &sol,int size);

    /*Compute no idle make span*/
    long int computeNIMS(std::vector<int> & sol);
    /*Compute no idle partial make span*/
    long int computeNIMS(std::vector<int> &sol, int size);
    /*Compute no idle make span iwthout computing the sums of machine 1 processing times*/
    long int computeNIMS(std::vector<int> &sol, long int nims);
    /*Compute no idle weighted tardiness*/
    long int computeNIWT(std::vector<int> &sol);
    long int computeNIWT(std::vector<int> &sol,int size);
    /*Compute no idle weighted earliness*/
    long int computeNIWE(std::vector<int> &sol);
    long int computeNIWE(std::vector<int> &sol,int size);
    /*Compute no idle earliness*/
    long int computeNIE(std::vector<int> &sol);
    long int computeNIE(std::vector<int> &sol,int size);
    /*Compute no idle tardiness*/
    long int computeNIT(std::vector<int> &sol);
    long int computeNIT(std::vector<int> &sol,int size);
    /*Compute no idle total completion time*/
    long int computeNITCT(std::vector< int > &sol);
    long int computeNITCT(std::vector< int > &sol,int size);
    long int computeNIWCT(std::vector< int > &sol);
    long int computeNIWCT(std::vector< int > &sol,int size);


    /*Compute Make Span with sequence depedent setup times*/
    long int computeSDSTMS(std::vector< int > &sol);
    long int computeSDSTMS(std::vector< int > &sol,int size);
    /*Compute sequence depedent setup times weighted earliness*/
    long int computeSDSTWT(std::vector<int> &sol);
    long int computeSDSTWT(std::vector<int> &sol,int size);
    /*Compute sequence depedent setup times weighted earliness*/
    long int computeSDSTWE(std::vector<int> &sol);
    long int computeSDSTWE(std::vector<int> &sol,int size);
    /*Compute sequence depedent setup times earliness*/
    long int computeSDSTE(std::vector<int> &sol);
    long int computeSDSTE(std::vector<int> &sol,int size);
    /*Compute sequence depedent setup times tardiness*/
    long int computeSDSTT(std::vector<int> &sol);
    long int computeSDSTT(std::vector<int> &sol,int size);
    /*Compute sequence depedent setup times total completion time*/
    long int computeSDSTTCT(std::vector< int > &sol);
    long int computeSDSTTCT(std::vector< int > &sol,int size);
    /*Compute sequence depedent setup times weighted completion time*/
    long int computeSDSTWCT(std::vector< int > &sol);
    long int computeSDSTWCT(std::vector< int > &sol,int size);


    /* Compute weighted tardines starting from an index*/
    long int computeWT(std::vector< int > & sol,std::vector<std::vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i);

    long int computeWT(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    //to compute WT and save the values for the job job
    void computeWTs(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    void setSilence(bool s);

    //Compute the earliest completion time starting from the beginning and starting from the end in order to apply taillard acceleration...
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size);
    void computeTails(std::vector<int> &sol, std::vector < std::vector< std::vector< int > > >& tails);

    void computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    int computeIdleTimeCoeff(std::vector<int>& prevJob, int job);

    const std::vector< std::vector < long int > > & getProcessingTimesMatrix() { return processingTimesMatrix; }



};

#endif