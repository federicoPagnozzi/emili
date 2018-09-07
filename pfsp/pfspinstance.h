//
//  Created by by Jérémie Dubois-Lacoste and Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.

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

    /**  Read write privates attributes : */
    int getNbJob();
    void setNbJob(int jobCount);

    int getNbMac();
    void setNbMac(int machienCount);

    /**  Allow the memory for the processing times matrix : */
    void allowMatrixMemory(int nbJ, int nbM);

    /**  Read\Write values in the matrix : */
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

    std::vector< std::vector < std::vector< int > > >& getSetUpTimes()
    {
        return setUpTimes;
    }

    /**  Read Data from a file : */
    bool readDataFromFile(char * fileName);

    /**  Read Data from file with the other format*/
    bool readDataFromFile(const std::string _filename);

    /**  Read Data from sequence dependent setup times file*/
    bool readSeqDepDataFromFile(char* filename);

    /** Compute weighted tardiness*/
    long int computeWT (std::vector< int > & sol);
    /** Compute partial weighted tardiness*/
    long int computeWT (std::vector< int > & sol, int size);
    long int computeWT(std::vector<int> &sol, std::vector<int>& makespans,int size);
    /**  Compute MakeSpan */
    long int computeMS (std::vector<int> & sol);
    /**  Compute partial MakeSpan*/
    long int computeMS(std::vector<int> &sol,int size);
    /**  Compute FlowTime*/
    long int computeFT(std::vector<int> & sol);
    long int computeFT(std::vector<int> &sol, int size);
    /** Compute weighted completion times*/
    /** Compute partial weighted completion times*/
    long int computeWCT (std::vector< int > & sol);    
    long int computeWCT (std::vector< int > & sol, int size);
    long int computeWCT(std::vector<int> &sol, std::vector<int>& makespans,int size);
    /** Compute weighted earliness*/
    /** Compute partial weighted earliness*/
    long int computeWE (std::vector< int > & sol);   
    long int computeWE (std::vector< int > & sol, int size);
    long int computeWE(std::vector<int> &sol, std::vector<int>& makespans,int size);
    /** Compute tardiness*/
    long int computeT(std::vector< int > & sol);
    /** Compute partial tardiness*/
    long int computeT(std::vector< int > & sol, int size);
    long int computeT(std::vector<int> &sol, std::vector<int>& makespans,int size);
    /** Compute earliness*/
    long int computeE (std::vector< int > & sol);
    /** Compute partial earliness*/
    long int computeE (std::vector< int > & sol, int size);
    long int computeE(std::vector<int> &sol, std::vector<int>& makespans,int size);
    /** Compute total completion time*/
    long int computeTCT(std::vector< int > &sol);
    long int computeTCT(std::vector< int > &sol,int size);
    long int computeTCT(std::vector<int> &sol, std::vector<int>& makespans,int size);

    /** Compute no wait make span*/
    long int computeNWMS(std::vector< int > & sol);
    /** Compute partial no wait make span*/
    long int computeNWMS(std::vector<int> & sol, int size);
    /** Compute no wait weighted tardiness*/
    long int computeNWWT(std::vector<int> &sol);
    long int computeNWWT(std::vector<int> &sol,int size);
    /** Compute no wait weighted earliness*/
    long int computeNWWE(std::vector<int> &sol);
    long int computeNWWE(std::vector<int> &sol,int size);
    /** Compute no wait earliness*/
    long int computeNWE(std::vector<int> &sol);
    long int computeNWE(std::vector<int> &sol,int size);
    /** Compute no wait tardiness*/
    long int computeNWT(std::vector<int> &sol);
    long int computeNWT(std::vector<int> &sol,int size);
    /** Compute no wait total completion time*/
    long int computeNWTCT(std::vector< int > &sol);
    long int computeNWTCT(std::vector< int > &sol,int size);
    long int computeNWWCT(std::vector< int > &sol);
    long int computeNWWCT(std::vector< int > &sol,int size);

    /** Compute no idle make span*/
    long int computeNIMS(std::vector<int> & sol);
    /** Compute no idle partial make span*/
    long int computeNIMS(std::vector<int> &sol, int size);
    /** Compute no idle make span without computing the sums of machine 1 processing times*/
    long int computeNIMS(std::vector<int> &sol, long int nims);
    /** Compute no idle weighted tardiness*/
    long int computeNIWT(std::vector<int> &sol);
    long int computeNIWT(std::vector<int> &sol,int size);
    /** Compute no idle weighted earliness*/
    long int computeNIWE(std::vector<int> &sol);
    long int computeNIWE(std::vector<int> &sol,int size);
    /** Compute no idle earliness*/
    long int computeNIE(std::vector<int> &sol);
    long int computeNIE(std::vector<int> &sol,int size);
    /** Compute no idle tardiness*/
    long int computeNIT(std::vector<int> &sol);
    long int computeNIT(std::vector<int> &sol,int size);
    /** Compute no idle total completion time*/
    long int computeNITCT(std::vector< int > &sol);
    long int computeNITCT(std::vector< int > &sol,int size);
    long int computeNIWCT(std::vector< int > &sol);
    long int computeNIWCT(std::vector< int > &sol,int size);


    /** Compute Make Span with sequence depedent setup times*/
    long int computeSDSTMS(std::vector< int > &sol);
    long int computeSDSTMS(std::vector< int > &sol,int size);
    /** Compute sequence depedent setup times weighted earliness*/
    long int computeSDSTWT(std::vector<int> &sol);
    long int computeSDSTWT(std::vector<int> &sol,int size);

    /** Compute sequence depedent setup times weighted earliness*/
    long int computeSDSTWE(std::vector<int> &sol);
    long int computeSDSTWE(std::vector<int> &sol,int size);

    /** Compute sequence depedent setup times earliness*/
    long int computeSDSTE(std::vector<int> &sol);
    long int computeSDSTE(std::vector<int> &sol,int size);

    /** Compute sequence depedent setup times tardiness*/
    long int computeSDSTT(std::vector<int> &sol);
    long int computeSDSTT(std::vector<int> &sol,int size);

    /** Compute sequence depedent setup times total completion time*/
    long int computeSDSTTCT(std::vector< int > &sol);
    long int computeSDSTTCT(std::vector< int > &sol,int size);

    /** Compute sequence depedent setup times weighted completion time*/
    long int computeSDSTWCT(std::vector< int > &sol);
    long int computeSDSTWCT(std::vector< int > &sol,int size);


    /** Compute Makespan LowerBound **/

    long int computeMSLB();

    /**  Compute weighted tardines starting from an index*/
    long int computeWT(std::vector< int > & sol,std::vector<std::vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i);
    //Compute weighted tardiness starting from the permutation and a partial calculation
    long int computeWT(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    //Starting from a solution computes the completion time for each machine ( prevjob) and the pre
    void computeWTs(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime);
    void setSilence(bool s);

    //Compute the earliest completion time starting from the beginning and starting from the end in order to apply taillard acceleration...
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    void computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size);
    void computeTails(std::vector<int> &sol, std::vector < std::vector< std::vector< int > > >& tails);
    void computeSDSTTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size);
    void computeSDSThead(std::vector<int> &sol,std::vector< std::vector < int > >& head, int size);
    void computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail);
    void computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size);
    int computeIdleTimeCoeff(std::vector<int>& prevJob, int job);
    void computeHead(std::vector<int>& sol,std::vector< std::vector< int > >& head, int njobs);
    //void updateHead(std::vector<int> &solution, int starting_point, std::vector < std::vector < int > >& head, std::vector<int>& makespans);

    const std::vector< std::vector < long int > > & getProcessingTimesMatrix() { return processingTimesMatrix; }   

};


#endif
