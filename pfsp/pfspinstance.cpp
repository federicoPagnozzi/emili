//
//  Created by by Jérémie Dubois-Lacoste and Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <exception>
#include <stdexcept>
#include "pfspinstance.h"
#include  <assert.h>

//#define ENABLE_SSE

#ifdef ENABLE_SSE

#ifdef __SSE__
#include <xmmintrin.h>
#include <emmintrin.h>
#endif
#endif

void show_head_and_tail(std::vector< std::vector < int > >& head,std::vector< std::vector < int > >& tail, int nbMac, int nbJob)
{
    for( int m = 0; m <= nbMac ; m++)
    {
        for(int n=0;n <= nbJob ; n++)
        {
            std::cout << " " << head[m][n];
        }
        std::cout << "\n";
    }

    std::cout << "\n";

    for( int m = 0; m <= nbMac ; m++)
    {
        for(int n=0;n <= nbJob ; n++)
        {
            std::cout << " " << tail[m][n];
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    std::cout << "\n";
}

PfspInstance::PfspInstance()
{
  silence = false;
}


PfspInstance::~PfspInstance()
{
}

int PfspInstance::getNbJob()
{
  return nbJob;
}

void PfspInstance::setNbJob(int jobCount)
{
    this->nbJob = jobCount;
}

int PfspInstance::getNbMac()
{
  return nbMac;
}

void PfspInstance::setNbMac(int machineCount)
{
    this->nbMac = machineCount;
}

void PfspInstance::setSilence(bool s)
{
    this->silence = s;
}

/**  Allow the memory for the processing times matrix : */
void PfspInstance::allowMatrixMemory(int nbJ, int nbM)
{
  processingTimesMatrix.resize(nbJ+1);          // 1st dimension

  for (int cpt = 0; cpt < nbJ+1; ++cpt)
    processingTimesMatrix[cpt].resize(nbM+1); // 2nd dimension

	dueDates.resize(nbJ+1);
	priority.resize(nbJ+1);
}


long int PfspInstance::getTime(int job, int machine)
{
  if (job == 0)
    return 0;
  else
  {
    if ((job < 1) || (job > nbJob) || (machine < 1) || (machine > nbMac))
      std::cout    << "ERROR. file:pfspInstance.cpp, method:getTime. Out of bound. job=" << job
          << ", machine=" << machine << std::endl;

    return processingTimesMatrix[job][machine];
  }
}


/**  Read the instance from file : */
bool PfspInstance::readDataFromFile(char * fileName)
{
	bool everythingOK = true;
	int j, m; // iterators
	long int readValue;
    std::string str;
    std::ifstream fileIn;

	char * aux2;
	char fileNameOK[100] = "";

	aux2 = (strrchr(fileName, '/'));

	if (aux2 == NULL)
		aux2 = fileName;
	else
		aux2 += 1;

	strcat(fileNameOK, aux2);
    if(!silence)
    {
    std::cout << "name : " << fileNameOK << std::endl;
    std::cout << "file : " << fileName << std::endl;
    }
	fileIn.open(fileName);

	if ( fileIn.is_open() ) {

        fileIn >> nbJob;
        fileIn >> nbMac;
        fileIn >> readValue;
        if(readValue == 12345)
        {
            std::string fname(fileName);
            return readDataFromFile(fname);
        }        
		allowMatrixMemory(nbJob, nbMac);
        if(!silence){
            std::cout << "File " << fileName << " is now open, start to read..." << std::endl;
            std::cout << "Number of jobs : " << nbJob << std::endl;
            std::cout << "Number of machines : " << nbMac << std::endl;
            std::cout << "Memory allowed." << std::endl;
            std::cout << "Start to read matrix..." << std::endl;
        }
		for (j = 1; j <= nbJob; ++j)
		{
            for (m = 1; m <= nbMac; ++m)
			{
                if(!(j==1 && m==1))
                {
				fileIn >> readValue; // The number of each machine, not important !
                }
                fileIn >> readValue; // Process Time

				processingTimesMatrix[j][m] = readValue;
			}
		}
        fileIn >> str; // this is not read

		for (j = 1; j <= nbJob; ++j)
		{
			fileIn >> readValue; // -1
			fileIn >> readValue;
			dueDates[j] = readValue;
			fileIn >> readValue; // -1
			fileIn >> readValue;
            priority[j] = readValue;
		}
        if(!silence)
        std::cout << "All is read from file." << std::endl;
		fileIn.close();
	}
	else
	{
        if(!silence)
        std::cout    << "ERROR. file:pfspInstance.cpp, method:readDataFromFile, "
                << "error while opening file " << fileName << std::endl;
        everythingOK = false;

	}

	return everythingOK;
}


bool PfspInstance::readSeqDepDataFromFile(char* fileName)
{
    bool everythingOK = true;
    int j, m; // iterators
    long int readValue;
    std::string str;
    std::ifstream fileIn;

    char * aux2;
    char fileNameOK[100] = "";

    aux2 = (strrchr(fileName, '/'));

    if (aux2 == NULL)
        aux2 = fileName;
    else
        aux2 += 1;

    strcat(fileNameOK, aux2);
    if(!silence)
    {
    std::cout << "name : " << fileNameOK << std::endl;
    std::cout << "file : " << fileName << std::endl;
    }
    fileIn.open(fileName);

    if ( fileIn.is_open() ) {

        fileIn >> nbJob;
        fileIn >> nbMac;
        fileIn >> readValue;
        if(readValue == 12345)
        {
            std::string fname(fileName);
            return readDataFromFile(fname);
        }
        allowMatrixMemory(nbJob, nbMac);

        if(!silence){
            std::cout << "File " << fileName << " is now open, start to read..." << std::endl;
            std::cout << "Number of jobs : " << nbJob << std::endl;
            std::cout << "Number of machines : " << nbMac << std::endl;
            std::cout << "Memory allowed." << std::endl;
            std::cout << "Start to read matrix..." << std::endl;
        }
        for (j = 1; j <= nbJob; ++j)
        {
            for (m = 1; m <= nbMac; ++m)
            {
                if(!(j==1 && m==1))
                {
                fileIn >> readValue; // The number of each machine, not important !
                }
                fileIn >> readValue; // Process Time

                processingTimesMatrix[j][m] = readValue;
            }
        }
        fileIn >> str; // this is not read
        if(str.compare("SSD")==0)
        {
            setUpTimes.resize(nbMac+1);
            for(int i=0;i<nbMac+1;i++){
                setUpTimes[i].resize(nbJob+1);
                for(int j=0;j<nbJob+1;j++)
                {
                    setUpTimes[i][j].resize(nbJob+1);
                }
            }

            for(int i=1; i< nbMac+1 ; i++)
            {
                fileIn >> str; // this is not read
                //std::cout << "M" << i << std::endl;
                for(int j=1;j<nbJob+1;j++)
                {
                    for(int k=1;k<nbJob+1;k++)
                    {
                        fileIn >> readValue; // -1
                        //fileIn >> readValue;
                        setUpTimes[i][j][k] = readValue;
                      //  std::cout << " " << readValue ;
                    }
                    //std::cout << "\n" ;
                }
            }



        }
        if(!fileIn.eof())
        {
        fileIn >> str; // this is not read
        if(str.compare("Reldue")==0)
        {
            for (j = 1; j <= nbJob; ++j)
           {
                fileIn >> readValue; // -1
                fileIn >> readValue;
                dueDates[j] = readValue;
                fileIn >> readValue; // -1
                fileIn >> readValue;
                priority[j] = readValue;
           }
        }
        }
        if(!silence)
        std::cout << "All is read from file." << std::endl;
        fileIn.close();
    }
    else
    {
        if(!silence)
        std::cout    << "ERROR. file:pfspInstance.cpp, method:readDataFromFile, "
                << "error while opening file " << fileName << std::endl;
        everythingOK = false;

    }

    return everythingOK;
}

bool PfspInstance::readDataFromFile(const std::string _fileName)
{
          std::string buffer;
          std::string::size_type start, end;
          std::ifstream inputFile(_fileName.data(), std::ios::in);
          // opening of the benchmark file
          //estd::cout << "file <- " << _fileName << std::endl;
          if (!inputFile.is_open())
              throw std::runtime_error("*** ERROR : Unable to open the benchmark file");
          // number of jobs (N)
          getline(inputFile, buffer, '\n');
          nbJob = atoi(buffer.data());
          // number of machines M
          getline(inputFile, buffer, '\n');
          nbMac = atoi(buffer.data());
          // initial and current seeds (not used)
          getline(inputFile, buffer, '\n');
          // processing times and due-dates
          // p = std::vector< std::vector<unsigned int> > (M,N);
          //p.resize(nbMac);
          allowMatrixMemory(nbJob, nbMac);
          /** for (unsigned int j=1 ; j<nbMac ; j++)
          {
              p[j].resize(N);
          }*/
          //d = std::vector<unsigned int> (nbJob);
          // for each job...
          for (unsigned int j=1 ; j<nbJob+1 ; j++)
          {
              // index of the job (<=> j)
              getline(inputFile, buffer, '\n');
              // due-date of the job j
              getline(inputFile, buffer, '\n');
              dueDates[j] = atoi(buffer.data());
              // processing times of the job j on each machine
              getline(inputFile, buffer, '\n');
              start = buffer.find_first_not_of(" ");
              for (unsigned int i=1 ; i<nbMac+1 ; i++)
              {
                  end = buffer.find_first_of(" ", start);
                  processingTimesMatrix[j][i] = atoi(buffer.substr(start, end-start).data());
                  start = buffer.find_first_not_of(" ", end);
              }
          }

          for (unsigned int j=1 ; j<nbJob+1 ; j++)
          {
            getline(inputFile, buffer, '\n');
            priority[j] = atoi(buffer.data());
          }

          // closing of the input file
          inputFile.close();
          if(!silence){
              std::cout << "File " << _fileName << " is now open, start to read..." << std::endl;
              std::cout << "Number of jobs : " << nbJob << std::endl;
              std::cout << "Number of machines : " << nbMac << std::endl;
              std::cout << "Memory allowed." << std::endl;
              std::cout << "Start to read matrix..." << std::endl;
          }
    return true;
}


#ifndef ENABLE_SSE
inline void computePartialMakespans( std::vector< int >& sol, std::vector< long int >& previousMachineEndTime,std::vector< std::vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    long int previousJobEndTime;
    int j, m;
    int jobNumber;
    /**  1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];
    }

    /**  others machines : */
    for ( m = 2; m <= nbMac; ++m )
    {
        previousMachineEndTime[1] +=
                processingTimesMatrix[sol[1]][m];
        previousJobEndTime = previousMachineEndTime[1];


        for ( j = 2; j <= nbJob; ++j )
        {
            jobNumber = sol[j];

            if ( previousMachineEndTime[j] > previousJobEndTime )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime += processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }
        }
    }
}

long int PfspInstance::computeMS(std::vector<int> &sol)
{
    int j, m;
    int jobNumber;

    /**  We need end times on previous machine : */
    std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /**  1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];
    }

    /**  others machines : */
    for ( m = 2; m <= nbMac; ++m )
    {
        previousMachineEndTime[1] += processingTimesMatrix[sol[1]][m];
        previousJobEndTime = previousMachineEndTime[1];


        for ( j = 2; j <= nbJob; ++j )
        {
            jobNumber = sol[j];

            if ( previousMachineEndTime[j] > previousJobEndTime )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime += processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }
        }
    }

    return previousMachineEndTime[nbJob];
}


long int PfspInstance::computeMS(std::vector<int> &sol,int size)
{
    int j, m;
    int jobNumber;

    // We need end times on previous machine :
    std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    // And the end time of the previous job, on the same machine :
    long int previousJobEndTime;
     //1st machine :
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= size; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];
    }

    // others machines :
    for ( m = 2; m <= nbMac; ++m )
    {
        previousMachineEndTime[1] += processingTimesMatrix[sol[1]][m];
        previousJobEndTime = previousMachineEndTime[1];


        for ( j = 2; j <= size; ++j )
        {
            jobNumber = sol[j];

            if ( previousMachineEndTime[j] > previousJobEndTime )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime += processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }
        }
    }

    return previousMachineEndTime[size];
}


#else



inline void computePartialMakespans(std::vector<int>& sol,std::vector< long >& previousMachineEndTime, std::vector<std::vector< long> >& pmat,int nbJob, int nbMac)
{

    /**  Permutation flowshop makespan computation using SSE instructions
     **/
    // Each sse register can contain 4 float so the computation is divided in groups of 4 jobs
    int r4 = nbJob%4;
    int lambda_number =  r4==0?nbJob/4:(nbJob/4+1); // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;
    if(r4>0)
    {
        for(int i=0;i<(4-r4);i++)
        {
            sol.push_back(0);
            previousMachineEndTime.push_back(0);
        }
    }
    int j=1;    
    std::vector<float> L(nbMac+1,0); // the makespan for each machine of the fourth job in the last group ( at the beginning is zero)
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);    
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 K = _mm_setzero_ps();

    for(int l = 0 ; l < lambda_number ; l++)
    {
        /**  At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations
        int j1 = sol[j],j2=sol[j+1],j3=sol[j+2],j4=sol[j+3];
        __m128 makespan, mc,mcw;
        makespan = _mm_set_ps(pmat[j4][1],pmat[j3][1],pmat[j2][1],pmat[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /** First machine
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1
        // Third add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a1,0,0,a0
        mc = _mm_and_ps(mc,mask0111);         // 0,0,0,a0
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1+a0
        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k

        K = _mm_shuffle_ps(makespan,makespan,0xFF);        

        /** The other machines
         *
         **/
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,pmat[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,pmat[j2][m],pmat[j1][m+1]);        // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)

        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,pmat[j3][m],
                           pmat[j2][m+1],pmat[j1][m+2]);        // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        //other rows
        for(m = 5; m <= nbMac ; m++)
        {
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(pmat[j4][m-3],pmat[j3][m-2],
                               pmat[j2][m-1],pmat[j1][m]);      // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            L[m-3] = res[3];
        }
                                                                //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        m = nbMac;
        mc = makespan;
        mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        mcw = _mm_set_ps(pmat[j4][m-2],
                           pmat[j3][m-1],pmat[j2][m],0);        // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,makespan);
        L[m-2] = res[3];
        // m - 2
        mc = makespan;
        mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        mcw = _mm_set_ps(pmat[j4][m-1],pmat[j3][m],0,0);        // mcw -> [ 0, 0 , T3,M,T4,M-1]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,makespan);
        L[m-1] = res[3];
        //m - 1
        mc = makespan;
        mc = _mm_and_ps(mc,mask0010);                           // mc -> [0 , 0, C3,M, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , 0, C3,M]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        mcw = _mm_set_ps(pmat[j4][m],0,0,0);                    // mcw -> [ 0, 0 , 0,T4,M]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,makespan);
        L[m] = res[3];
        previousMachineEndTime[j] = res[0];
        previousMachineEndTime[j+1] = res[1];
        previousMachineEndTime[j+2] = res[2];
        previousMachineEndTime[j+3] = res[3];
        j+=4;

        /**  At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
    if(r4>0)
    {
        for(int i=0;i<(4-r4);i++)
        {
            sol.pop_back();
            previousMachineEndTime.pop_back();
        }
    }

}

long int PfspInstance::computeMS(std::vector<int> &sol)
{
    return computeMS(sol,nbJob);
}

long int PfspInstance::computeMS(std::vector<int> &sol, int size)
{
    /**  Permutation flowshop makespan computation using SSE instructions
     **/
    // Each sse register can contain 4 float so the computation is divided in groups of 4 jobs

    int r4 = size%4;
    int lambda_number =  r4==0?size/4:(size/4+1); // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;
    if(r4>0)
    {
        for(int i=0;i<(4-r4);i++)
        {
            sol.push_back(0);
        }
    }
    int j=1;
    std::vector<float> L(nbMac+1,0); // the makespan for each machine of the fourth job in the last group ( at the beginning is zero)
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 K = _mm_setzero_ps();

    for(int l = 0 ; l < lambda_number ; l++)
    {
        /**  At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations
        int j1 = sol[j],j2=sol[j+1],j3=sol[j+2],j4=sol[j+3];
        __m128 makespan, mc,mcw;
        makespan = _mm_set_ps(processingTimesMatrix[j4][1],
                processingTimesMatrix[j3][1],
                processingTimesMatrix[j2][1],
                processingTimesMatrix[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /** First machine
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1
        // Third add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a1,0,0,a0
        mc = _mm_and_ps(mc,mask0111);         // 0,0,0,a0
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1+a0
        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k

        K = _mm_shuffle_ps(makespan,makespan,0xFF);

        /** The other machines
         *
         **/
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,processingTimesMatrix[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,processingTimesMatrix[j2][m],
                         processingTimesMatrix[j1][m+1]);       // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)

        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,processingTimesMatrix[j3][m],
                           processingTimesMatrix[j2][m+1],
                processingTimesMatrix[j1][m+2]);                // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        //other rows
        for(m = 5; m <= nbMac ; m++)
        {
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(processingTimesMatrix[j4][m-3],
                    processingTimesMatrix[j3][m-2],
                               processingTimesMatrix[j2][m-1],
                    processingTimesMatrix[j1][m]);              // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            L[m-3] = res[3];
        }
                                                                //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        m = nbMac;
        mc = makespan;
        mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-2],
                           processingTimesMatrix[j3][m-1],
                processingTimesMatrix[j2][m],0);                // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,makespan);
        L[m-2] = res[3];
        // m - 2
        mc = makespan;
        mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-1],
                processingTimesMatrix[j3][m],0,0);              // mcw -> [ 0, 0 , T3,M,T4,M-1]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,makespan);
        L[m-1] = res[3];
        //m - 1
        mc = makespan;
        mc = _mm_and_ps(mc,mask0010);                           // mc -> [0 , 0, C3,M, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , 0, C3,M]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m],
                         0,0,0);                                // mcw -> [ 0, 0 , 0,T4,M]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,makespan);
        L[m] = res[3];
        j+=4;

        /**  At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
    long int result = res[3];
    if(r4>0)
    {
        for(int i=0;i<(4-r4);i++)
        {
            sol.pop_back();

        }
        result = res[r4-1];
    }
    return result;
}

inline void computeHEADandTAIL(std::vector<int> &sol,std::vector< std::vector < int > >& head,std::vector< std::vector < int > >& tail,std::vector<std::vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    /**  Permutation flowshop Tail and Head matrices computation using SSE instructions
     **/
    // Each sse register can contain 4 float so the computation is divided in groups of 4 jobs

    int r4 = nbJob%4;
    int lambda_number =  nbJob/4; // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;

    int j=1;
    int tj = nbJob;
    std::vector<float> L(nbMac+1,0); // the makespan for each machine of the fourth job in the last group ( at the beginning is zero)
    std::vector<float> tL(nbMac+1,0);
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 K = _mm_setzero_ps();
    __m128 tK = _mm_setzero_ps();
    for(int l = 0 ; l < lambda_number ; l++)
    {
        /**  At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations
        int j1 = sol[j],j2=sol[j+1],j3=sol[j+2],j4=sol[j+3];
        int tj1 = sol[tj],tj2=sol[tj-1],tj3=sol[tj-2],tj4=sol[tj-3];
        __m128 makespan, mc,mcw;
        __m128 tailspan,tc,tcw;
        makespan = _mm_set_ps(processingTimesMatrix[j4][1],
                processingTimesMatrix[j3][1],
                processingTimesMatrix[j2][1],
                processingTimesMatrix[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /** First machine HEAD
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1
        // Third add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a1,0,0,a0
        mc = _mm_and_ps(mc,mask0111);         // 0,0,0,a0
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1+a0
        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k

        K = _mm_shuffle_ps(makespan,makespan,0xFF);
        _mm_store_ps(res,makespan);

        head[1][j] = res[0];
        head[1][j+1] = res[1];
        head[1][j+2] = res[2];
        head[1][j+3] = res[3];
        /** First machine TAIL
         *
         * */

        tailspan = _mm_set_ps(processingTimesMatrix[tj4][1],
                processingTimesMatrix[tj3][1],
                processingTimesMatrix[tj2][1],
                processingTimesMatrix[tj1][1]); // load the values in the registers
        tc = tailspan; // copy the value in another register
                                              // a3,a2,a1,a0
        // first add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a0,a3,a2,a1
        tc = _mm_and_ps(tc,mask0111);         // 0,a3,a2,a1
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1,a0+a1
        // Second add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a1,0,a3,a2
        tc = _mm_and_ps(tc,mask0111);         // 0,0,a3,a2
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1+a3,a2+a0+a1
        // Third add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a2,0,0,a3
        tc = _mm_and_ps(tc,mask0111);         // 0,0,0,a3
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a1+a2+a3,a1+a0+a2+a3
        tailspan = _mm_add_ps(tailspan,tK);    // a3+k,a3+a2+k,a1+a2+a3+k,a1+a0+a2+a3+k

        tK = _mm_shuffle_ps(tailspan,tailspan,0xFF);
        _mm_store_ps(res,tailspan);

        tail[nbMac][tj] = res[0];
        tail[nbMac][tj-1] = res[1];
        tail[nbMac][tj-2] = res[2];
        tail[nbMac][tj-3] = res[3];

        /** The other machines
         *
         **/

        /**  Filling the pipeline HEAD
         *
         * */
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,processingTimesMatrix[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,makespan);
        head[2][j] = res[0];

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,processingTimesMatrix[j2][m],
                         processingTimesMatrix[j1][m+1]);       // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,makespan);
        head[3][j] = res[0];
        head[2][j+1] = res[1];
        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,processingTimesMatrix[j3][m],
                           processingTimesMatrix[j2][m+1],
                processingTimesMatrix[j1][m+2]);                // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,makespan);
        head[4][j] = res[0];
        head[3][j+1] = res[1];
        head[2][j+2] = res[2];
        /**  Filling the pipeline TAIL
         *
         * */
        int tm = nbMac-1;
                                                                //tailspan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        tcw = _mm_set_ps(0,0,0,tL[tm]);                           // tcw -> [L2 ,0,0,0]
        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,0,processingTimesMatrix[tj1][tm]);  // tcw -> [T1,2,0,0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm][tj] = res[0];

        // second row
        tcw = _mm_set_ps(0,0,0,tL[tm-1]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1000);                           // tc -> [C1,2 , 0   , 0 , 0]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0   , C1,2, 0 , 0]
        tcw = _mm_add_ps(tcw,tc);                               // tcw ->[L3   , C1,2, 0 , 0]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,processingTimesMatrix[tj2][tm],
                         processingTimesMatrix[tj1][tm-1]);       // tcw -> [ T1,3, T2,2 , 0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,tailspan);
        tail[tm-1][tj] = res[0];
        tail[tm][tj-1] = res[1];
        // third row
        tcw = _mm_set_ps(0,0,0,tL[tm-2]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1100);                           // tc -> [C1,3 , C2,2, 0, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,3, C2,2, 0 ]
        tcw = _mm_add_ps(tcw,tc);                               // tcw -> [L4, C1,3, C2,2  , 0 ]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        tcw = _mm_set_ps(0,processingTimesMatrix[tj3][tm],
                           processingTimesMatrix[tj2][tm-1],
                processingTimesMatrix[tj1][tm-2]);                // tcw -> [ T1,4, T2,3 , T3,2,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj] = res[0];
        tail[tm-1][tj-1] = res[1];
        tail[tm][tj-2] = res[2];
        //other rows
        tm = tm-3;
        for(m = 5; m <= nbMac ; m++)
        {
            /** HEAD
             *
             * */
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(processingTimesMatrix[j4][m-3],
                    processingTimesMatrix[j3][m-2],
                               processingTimesMatrix[j2][m-1],
                    processingTimesMatrix[j1][m]);              // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            head[m][j] = res[0];
            head[m-1][j+1] = res[1];
            head[m-2][j+2] = res[2];
            head[m-3][j+3] = res[3];
            L[m-3] = res[3];

            /** TAIL
             *
             * */

            // m row
            tcw = _mm_set_ps(0,0,0,tL[tm]);                       // setup vec for compares
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                       // tc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                    // tc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            tcw = _mm_add_ps(tcw,tc);                           // tcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            tailspan = _mm_max_ps(tcw,tailspan);                // tailspan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            tcw = _mm_set_ps(processingTimesMatrix[tj4][tm+3],
                    processingTimesMatrix[tj3][tm+2],
                               processingTimesMatrix[tj2][tm+1],
                    processingTimesMatrix[tj1][tm]);              // tcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            tailspan = _mm_add_ps(tailspan,tcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,tailspan);
            tail[tm][tj] = res[0];
            tail[tm+1][tj-1] = res[1];
            tail[tm+2][tj-2] = res[2];
            tail[tm+3][tj-3] = res[3];
            tL[tm] = res[3];
            tm--;
        }
        /**
         * HEAD finishing...
         * */
                                                        //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        m = nbMac;
        mc = makespan;
        mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-2],
                           processingTimesMatrix[j3][m-1],
                processingTimesMatrix[j2][m],0);                // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,makespan);
        head[m][j+1] = res[1];
        head[m-1][j+2] = res[2];
        head[m-2][j+3] = res[3];
        L[m-2] = res[3];
        // m - 2
        mc = makespan;
        mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m-1],
                processingTimesMatrix[j3][m],0,0);              // mcw -> [ 0, 0 , T3,M,T4,M-1]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,makespan);
        head[m][j+2] = res[2];
        head[m-1][j+3] = res[3];
        L[m-1] = res[3];
        //m - 1
        mc = makespan;
        mc = _mm_and_ps(mc,mask0010);                           // mc -> [0 , 0, C3,M, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , 0, C3,M]
        makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        mcw = _mm_set_ps(processingTimesMatrix[j4][m],
                         0,0,0);                                // mcw -> [ 0, 0 , 0,T4,M]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,makespan);
        head[m][j+3] = res[3];
        L[m] = res[3];
        j+=4;

        /** TAIL finishing
         *
         * */
        // m - 3
        tm = 3;
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1110);                           // tc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,M , C2,M-1, C3,M-2]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm],
                           processingTimesMatrix[tj3][tm-1],
                processingTimesMatrix[tj2][tm-2],0);                // tcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //tailspan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-1] = res[1];
        tail[tm-1][tj-2] = res[2];
        tail[tm][tj-3] = res[3];
        tL[tm] = res[3];
        // m - 2
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0110);                           // tc -> [0 , C2,M, C3,M-1, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , C2,M, C3,M-1]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-1],
                processingTimesMatrix[tj3][tm-2],0,0);              // tcw -> [ 0, 0 , T3,M,T4,M-1]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-2] = res[2];
        tail[tm-1][tj-3] = res[3];
        tL[tm-1] = res[3];
        //m - 1
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0010);                           // tc -> [0 , 0, C3,M, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , 0, C3,M]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-2],
                         0,0,0);                                // tcw -> [ 0, 0 , 0,T4,M]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-3] = res[3];
        tL[tm] = res[3];
        tj-=4;
        /**  At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
// What if the number of jobs cannot be divide by 4?
    if(r4 > 0)
    {
        //Initializations
        int j1 = sol[j],j2= r4>=2?sol[j+1]:0,j3=r4>=3?sol[j+2]:0;
        int tj1 = sol[tj],tj2=r4>=2?sol[tj-1]:0,tj3=r4>=3?sol[tj-2]:0;
        __m128 makespan, mc,mcw;
        __m128 tailspan,tc,tcw;
        makespan = _mm_set_ps(0,
                processingTimesMatrix[j3][1],
                processingTimesMatrix[j2][1],
                processingTimesMatrix[j1][1]); // load the values in the registers
        mc = makespan; // copy the value in another register

        /** First machine HEAD
         *
         * */
        // first add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a3,a0,a1,a2
        mc = _mm_and_ps(mc,mask0111);         // 0,a0,a1,a2
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1,a2+a3
        // Second add
        mc = _mm_shuffle_ps(mc,mc,0x93);      // a2,0,a0,a1
        mc = _mm_and_ps(mc,mask0111);         // 0,0,a0,a1
        makespan = _mm_add_ps(makespan,mc);   // a0,a0+a1,a2+a1+a0,a2+a3+a1

        makespan = _mm_add_ps(makespan,K);    // a0+k,a0+a1+k,a2+a1+a0+k,a2+a3+a1+a0+k


        _mm_store_ps(res,makespan);

        head[1][j] = res[0];
        if(r4 > 1)
        {
            head[1][j+1] = res[1];
            if(r4 > 2)
            head[1][j+2] = res[2];
        }
        /** First machine TAIL
         *
         * */

        tailspan = _mm_set_ps(0,
                processingTimesMatrix[tj3][1],
                processingTimesMatrix[tj2][1],
                processingTimesMatrix[tj1][1]); // load the values in the registers
        tc = tailspan; // copy the value in another register
                                              // a3,a2,a1,a0
        // first add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a0,a3,a2,a1
        tc = _mm_and_ps(tc,mask0111);         // 0,a3,a2,a1
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1,a0+a1
        // Second add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a1,0,a3,a2
        tc = _mm_and_ps(tc,mask0111);         // 0,0,a3,a2
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1+a3,a2+a0+a1

        tailspan = _mm_add_ps(tailspan,tK);    // a3+k,a3+a2+k,a1+a2+a3+k,a1+a0+a2+a3+k

        _mm_store_ps(res,tailspan);

        tail[nbMac][tj] = res[0];
        if(r4 > 1)
        {
            tail[nbMac][tj-1] = res[1];
            if(r4 > 2)
            tail[nbMac][tj-2] = res[2];
        }

        /** The other machines
         *
         **/

        /**  Filling the pipeline HEAD
         *
         * */
        int m=2;
                                                                //makespan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        mcw = _mm_set_ps(0,0,0,L[2]);                           // mcw -> [L2 ,0,0,0]
        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,0,processingTimesMatrix[j1][m]);                    // mcw -> [T1,2,0,0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,makespan);
        head[2][j] = res[0];

        // second row
        mcw = _mm_set_ps(0,0,0,L[3]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1000);                           // mc -> [C1,2 , 0   , 0 , 0]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0   , C1,2, 0 , 0]
        mcw = _mm_add_ps(mcw,mc);                               // mcw ->[L3   , C1,2, 0 , 0]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        mcw = _mm_set_ps(0,0,processingTimesMatrix[j2][m],
                         processingTimesMatrix[j1][m+1]);       // mcw -> [ T1,3, T2,2 , 0,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,makespan);
        head[3][j] = res[0];
        if(r4>1)
        head[2][j+1] = res[1];
        // third row
        mcw = _mm_set_ps(0,0,0,L[4]);                           // setup vec for compares
        mc = makespan;
        mc = _mm_and_ps(mc,mask1100);                           // mc -> [C1,3 , C2,2, 0, 0 ]
        mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,3, C2,2, 0 ]
        mcw = _mm_add_ps(mcw,mc);                               // mcw -> [L4, C1,3, C2,2  , 0 ]

        makespan = _mm_max_ps(mcw,makespan);                    // makespan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        mcw = _mm_set_ps(0,processingTimesMatrix[j3][m],
                           processingTimesMatrix[j2][m+1],
                processingTimesMatrix[j1][m+2]);                // mcw -> [ T1,4, T2,3 , T3,2,0]
        makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,makespan);
        head[4][j] = res[0];
        if(r4>1)
        head[3][j+1] = res[1];
        if(r4>2)
        head[2][j+2] = res[2];
        /**  Filling the pipeline TAIL
         *
         * */
        int tm = nbMac-1;
                                                                //tailspan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        tcw = _mm_set_ps(0,0,0,tL[tm]);                           // tcw -> [L2 ,0,0,0]
        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,0,processingTimesMatrix[tj1][tm]);  // tcw -> [T1,2,0,0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm][tj] = res[0];

        // second row
        tcw = _mm_set_ps(0,0,0,tL[tm-1]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1000);                           // tc -> [C1,2 , 0   , 0 , 0]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0   , C1,2, 0 , 0]
        tcw = _mm_add_ps(tcw,tc);                               // tcw ->[L3   , C1,2, 0 , 0]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,processingTimesMatrix[tj2][tm],
                         processingTimesMatrix[tj1][tm-1]);       // tcw -> [ T1,3, T2,2 , 0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,tailspan);
        tail[tm-1][tj] = res[0];
        if(r4>1)
        tail[tm][tj-1] = res[1];
        // third row
        tcw = _mm_set_ps(0,0,0,tL[tm-2]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1100);                           // tc -> [C1,3 , C2,2, 0, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,3, C2,2, 0 ]
        tcw = _mm_add_ps(tcw,tc);                               // tcw -> [L4, C1,3, C2,2  , 0 ]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        tcw = _mm_set_ps(0,processingTimesMatrix[tj3][tm],
                           processingTimesMatrix[tj2][tm-1],
                processingTimesMatrix[tj1][tm-2]);                // tcw -> [ T1,4, T2,3 , T3,2,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj] = res[0];
        if(r4>1)
        tail[tm-1][tj-1] = res[1];
        if(r4>2)
        tail[tm][tj-2] = res[2];
        //other rows
        tm = tm-3;
        for(m = 5; m <= nbMac ; m++)
        {
            /** HEAD
             *
             * */
            // m row
            mcw = _mm_set_ps(0,0,0,L[m]);                       // setup vec for compares
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                       // mc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                    // mc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            mcw = _mm_add_ps(mcw,mc);                           // mcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            makespan = _mm_max_ps(mcw,makespan);                // makespan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            mcw = _mm_set_ps(0,
                    processingTimesMatrix[j3][m-2],
                               processingTimesMatrix[j2][m-1],
                    processingTimesMatrix[j1][m]);              // mcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            makespan = _mm_add_ps(makespan,mcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //makespan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,makespan);
            head[m][j] = res[0];
            if(r4>1)
            head[m-1][j+1] = res[1];
            if(r4>2)
            head[m-2][j+2] = res[2];
            /** TAIL
             *
             * */
            tcw = _mm_set_ps(0,0,0,tL[tm]);                       // setup vec for compares
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                       // tc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                    // tc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            tcw = _mm_add_ps(tcw,tc);                           // tcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            tailspan = _mm_max_ps(tcw,tailspan);                // tailspan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            tcw = _mm_set_ps(0,
                    processingTimesMatrix[tj3][tm+2],
                               processingTimesMatrix[tj2][tm+1],
                    processingTimesMatrix[tj1][tm]);              // tcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            tailspan = _mm_add_ps(tailspan,tcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,tailspan);
            tail[tm][tj] = res[0];
            if(r4>1)
            tail[tm+1][tj-1] = res[1];
            if(r4>2)
            tail[tm+2][tj-2] = res[2];
            tm--;
        }
        /**
         * HEAD finishing...
         * */
                                                        //makespan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        if(r4>1)
        {
            m = nbMac;
            mc = makespan;
            mc = _mm_and_ps(mc,mask1110);                           // mc -> [C1,M , C2,M-1, C3,M-2, 0 ]
            mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, C1,M , C2,M-1, C3,M-2]
            makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
            mcw = _mm_set_ps(0,
                             processingTimesMatrix[j3][m-1],
                    processingTimesMatrix[j2][m],0);                // mcw -> [ 0, T2,M , T3,M-1,T4,M-2]
            makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
            //makespan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
            _mm_store_ps(res,makespan);
            head[m][j+1] = res[1];
            if(r4>2)
            {
                head[m-1][j+2] = res[2];
                // m - 2
                mc = makespan;
                mc = _mm_and_ps(mc,mask0110);                           // mc -> [0 , C2,M, C3,M-1, 0 ]
                mc = _mm_shuffle_ps(mc,mc,0x93);                        // mc -> [ 0, 0 , C2,M, C3,M-1]
                makespan = _mm_max_ps(mc,makespan);                     // makespan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
                mcw = _mm_set_ps(0,
                        processingTimesMatrix[j3][m],0,0);              // mcw -> [ 0, 0 , T3,M,T4,M-1]
                makespan = _mm_add_ps(makespan,mcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                //makespan -> [ C1,M, C2,M, C3,M, C4,M-1]
                _mm_store_ps(res,makespan);
                head[m][j+2] = res[2];
            }
        }
        /** TAIL finishing
         *
         * */
        // m - 3
        if(r4 > 1)
        {
            tm = 3;
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                           // tc -> [C1,M , C2,M-1, C3,M-2, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,M , C2,M-1, C3,M-2]
            tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
            tcw = _mm_set_ps(0,
                             processingTimesMatrix[tj3][tm-1],
                    processingTimesMatrix[tj2][tm-2],0);                // tcw -> [ 0, T2,M , T3,M-1,T4,M-2]
            tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
            //tailspan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
            _mm_store_ps(res,tailspan);
            tail[tm-2][tj-1] = res[1];
            if(r4>2)
            {
                tail[tm-1][tj-2] = res[2];
                // m - 2
                tc = tailspan;
                tc = _mm_and_ps(tc,mask0110);                           // tc -> [0 , C2,M, C3,M-1, 0 ]
                tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , C2,M, C3,M-1]
                tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
                tcw = _mm_set_ps(0,
                        processingTimesMatrix[tj3][tm-2],0,0);              // tcw -> [ 0, 0 , T3,M,T4,M-1]
                tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                //tailspan -> [ C1,M, C2,M, C3,M, C4,M-1]
                _mm_store_ps(res,tailspan);
                tail[tm-2][tj-2] = res[2];
            }
        }
        /**  At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }

}

inline void computeTAIL(std::vector<int> &sol,std::vector< std::vector < int > >& tail,std::vector<std::vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    /**  Permutation flowshop Tail and Head matrices computation using SSE instructions
     **/
    // Each sse register can contain 4 float so the computation is divided in groups of 4 jobs

    int r4 = nbJob%4;
    int lambda_number =  r4==0?nbJob/4:(nbJob/4+1); // if ( nbjob%4==0) lambda_number = nbjob/4 else lambda_number = nbjob/4+1;
    if(r4>0)
    {
        for(int i=0;i<r4;i++)
        {
            sol.push_back(0);
        }
    }

    int tj = nbJob;
// the tailspan for each machine of the fourth job in the last group ( at the beginning is zero)
    std::vector<float> tL(nbMac+1,0);
    float res[4] __attribute__((aligned(16)));
    int* k = (int*)res;
        k[0] = 0xffffffff;
        k[1] = 0xffffffff;
        k[2] = 0xffffffff;
        k[3] = 0;
    __m128 mask1110 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0xffffffff;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1100 = _mm_load_ps(res);
    k[0] = 0xffffffff;
    k[1] = 0;
    k[2] = 0;
    k[3] = 0;
    __m128 mask1000 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0xffffffff;
    __m128 mask0111 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0xffffffff;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0110 = _mm_load_ps(res);
    k[0] = 0;
    k[1] = 0;
    k[2] = 0xffffffff;
    k[3] = 0;
    __m128 mask0010 = _mm_load_ps(res);
    __m128 tK = _mm_setzero_ps();
    for(int l = 0 ; l < lambda_number ; l++)
    {
        /**  At the beginning
         * J1   J2   J3   J4
        K  T1,1 T2,1 T3,1 T4,1
        L2 T1,2 T2,2 T3,2 T4,2
        L3 T1,3 T2,3 T3,3 T4,3
        .. ..   ..   ..   ..
        LM T1,M T2,M T3,M T4,M

        res[ 0 , 0 , 0 , 0]

        */
        //Initializations

        int tj1 = sol[tj],tj2=sol[tj-1],tj3=sol[tj-2],tj4=sol[tj-3];
        __m128 tailspan,tc,tcw;
        /** First machine TAIL
         *
         * */

        tailspan = _mm_set_ps(processingTimesMatrix[tj4][1],
                processingTimesMatrix[tj3][1],
                processingTimesMatrix[tj2][1],
                processingTimesMatrix[tj1][1]); // load the values in the registers
        tc = tailspan; // copy the value in another register
                                              // a3,a2,a1,a0
        // first add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a0,a3,a2,a1
        tc = _mm_and_ps(tc,mask0111);         // 0,a3,a2,a1
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1,a0+a1
        // Second add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a1,0,a3,a2
        tc = _mm_and_ps(tc,mask0111);         // 0,0,a3,a2
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a2+a1+a3,a2+a0+a1
        // Third add
        tc = _mm_shuffle_ps(tc,tc,0x93);      // a2,0,0,a3
        tc = _mm_and_ps(tc,mask0111);         // 0,0,0,a3
        tailspan = _mm_add_ps(tailspan,tc);   // a3,a3+a2,a1+a2+a3,a1+a0+a2+a3
        tailspan = _mm_add_ps(tailspan,tK);    // a3+k,a3+a2+k,a1+a2+a3+k,a1+a0+a2+a3+k

        tK = _mm_shuffle_ps(tailspan,tailspan,0xFF);
        _mm_store_ps(res,tailspan);

        tail[nbMac][tj] = res[0];
        tail[nbMac][tj-1] = res[1];
        tail[nbMac][tj-2] = res[2];
        tail[nbMac][tj-3] = res[3];



        /** The other machines
         *
         **/
        /**  Filling the pipeline TAIL
         *
         * */
        int tm = nbMac-1;
                                                                //tailspan -> [ C1,1, C2,1 , C3,1 , C4,1]
        // first row
        tcw = _mm_set_ps(0,0,0,tL[tm]);                           // tcw -> [L2 ,0,0,0]
        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L2,C1,1),C2,1 , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,0,processingTimesMatrix[tj1][tm]);  // tcw -> [T1,2,0,0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,2, C2,1 , C3,1 , C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm][tj] = res[0];

        // second row
        tcw = _mm_set_ps(0,0,0,tL[tm-1]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1000);                           // tc -> [C1,2 , 0   , 0 , 0]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0   , C1,2, 0 , 0]
        tcw = _mm_add_ps(tcw,tc);                               // tcw ->[L3   , C1,2, 0 , 0]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L3,C1,2),max(C1,2 , C2,1) , C3,1 , C4,1]
        tcw = _mm_set_ps(0,0,processingTimesMatrix[tj2][tm],
                         processingTimesMatrix[tj1][tm-1]);       // tcw -> [ T1,3, T2,2 , 0,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
        _mm_store_ps(res,tailspan);
        tail[tm-1][tj] = res[0];
        tail[tm][tj-1] = res[1];
        // third row
        tcw = _mm_set_ps(0,0,0,tL[tm-2]);                           // setup vec for compares
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1100);                           // tc -> [C1,3 , C2,2, 0, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,3, C2,2, 0 ]
        tcw = _mm_add_ps(tcw,tc);                               // tcw -> [L4, C1,3, C2,2  , 0 ]

        tailspan = _mm_max_ps(tcw,tailspan);                    // tailspan -> [ max(L4,C1,2),max(C1,3 , C2,2) , max( C2,2, C3,1 ) , C4,1]
        tcw = _mm_set_ps(0,processingTimesMatrix[tj3][tm],
                           processingTimesMatrix[tj2][tm-1],
                processingTimesMatrix[tj1][tm-2]);                // tcw -> [ T1,4, T2,3 , T3,2,0]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,4, C2,3, C3,2, C4,1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj] = res[0];
        tail[tm-1][tj-1] = res[1];
        tail[tm][tj-2] = res[2];
        //other rows
        for(int m = tm-3; m >= 1 ; m--)
        {
            // m row
            tcw = _mm_set_ps(0,0,0,tL[m]);                       // setup vec for compares
            tc = tailspan;
            tc = _mm_and_ps(tc,mask1110);                       // tc -> [C1,m-1 , C2,m-2, C3,m-3, 0 ]
            tc = _mm_shuffle_ps(tc,tc,0x93);                    // tc -> [ 0, C1,m-1 , C2,m-2, C3,m-3]
            tcw = _mm_add_ps(tcw,tc);                           // tcw ->[Lm, C1,m-1 , C2,m-2, C3,m-3]

            tailspan = _mm_max_ps(tcw,tailspan);                // tailspan -> [ max(Lm,C1,m-1),max(C1,m-1 , C2,m-2) , max( C2,m-2, C3,m-3 ) , max(C3,m-3, C4,m-4) ]
            tcw = _mm_set_ps(processingTimesMatrix[tj4][m+3],
                    processingTimesMatrix[tj3][m+2],
                               processingTimesMatrix[tj2][m+1],
                    processingTimesMatrix[tj1][m]);              // tcw -> [ T1,m, T2,m-1 , T3,m-2,T4,m-3]
            tailspan = _mm_add_ps(tailspan,tcw);                // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,m, C2,m-1, C3,m-2, C4,m-3]
            _mm_store_ps(res,tailspan);
            tail[m][tj] = res[0];
            tail[m+1][tj-1] = res[1];
            tail[m+2][tj-2] = res[2];
            tail[m+3][tj-3] = res[3];
            tL[m] = res[3];
        }
                                                                //tailspan -> [C1,M , C2,M-1, C3,M-2, C4,M-3]
        // m - 3
        tm = 3;
        tc = tailspan;
        tc = _mm_and_ps(tc,mask1110);                           // tc -> [C1,M , C2,M-1, C3,M-2, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, C1,M , C2,M-1, C3,M-2]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M,max(C1,M , C2,M-1) , max( C2,M-1, C3,M-2 ) , max(C3,M-2, C4,M-3) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm],
                           processingTimesMatrix[tj3][tm-1],
                processingTimesMatrix[tj2][tm-2],0);                // tcw -> [ 0, T2,M , T3,M-1,T4,M-2]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,M-1 , Cj-1,M)
                                                                //tailspan -> [ C1,M, C2,M, C3,M-1, C4,M-2]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-1] = res[1];
        tail[tm-1][tj-2] = res[2];
        tail[tm][tj-3] = res[3];
        tL[tm] = res[3];
        // m - 2
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0110);                           // tc -> [0 , C2,M, C3,M-1, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , C2,M, C3,M-1]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , max( C2,M, C3,M-1 ) , max(C3,M-1, C4,M-2) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-1],
                processingTimesMatrix[tj3][tm-2],0,0);              // tcw -> [ 0, 0 , T3,M,T4,M-1]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M-1]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-2] = res[2];
        tail[tm-1][tj-3] = res[3];
        tL[tm-1] = res[3];
        //m - 1
        tc = tailspan;
        tc = _mm_and_ps(tc,mask0010);                           // tc -> [0 , 0, C3,M, 0 ]
        tc = _mm_shuffle_ps(tc,tc,0x93);                        // tc -> [ 0, 0 , 0, C3,M]
        tailspan = _mm_max_ps(tc,tailspan);                     // tailspan -> [ C1,M ,C2,M , C3,M , max(C3,M, C4,M-1) ]
        tcw = _mm_set_ps(processingTimesMatrix[tj4][tm-2],
                         0,0,0);                                // tcw -> [ 0, 0 , 0,T4,M]
        tailspan = _mm_add_ps(tailspan,tcw);                    // Tjm + max(Cj,m-1 , Cj-1,m)
                                                                //tailspan -> [ C1,M, C2,M, C3,M, C4,M]
        _mm_store_ps(res,tailspan);
        tail[tm-2][tj-3] = res[3];
        tL[tm] = res[3];
        tj-=4;

        /**  At the end
         * J1   J2   J3   J4
           T1,1 T2,1 T3,1 T4,1  K
           T1,2 T2,2 T3,2 T4,2  L2
           T1,3 T2,3 T3,3 T4,3  L3
           ..   ..   ..   ..    ..
           T1,M T2,M T3,M T4,M  LM

        res[ C1,M , C2,M , C3,M , C4,M]

        */
    }
    long int result = res[3];
    if(r4>0)
    {
        for(int i=0;i<r4;i++)
        {
            sol.pop_back();

        }
        result = res[r4-1];
    }
}


#endif
/**  Compute the weighted tardiness of a given solution */

long int PfspInstance::computeWT(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
    std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    /** !!! It could be implemented using SIMD instructions... */
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}
/** compute partial weighted tardiness*/
long int PfspInstance::computeWT(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
    std::vector< long int > previousMachineEndTime ( nbJob + 1 ,0);
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);


    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);
    }

    return wt;

}



long int PfspInstance::getDueDate(int job)
{
    return dueDates[job];
}

long int PfspInstance::getPriority(int job)
{
    return priority[job];
}


/**  Compute the weighted tardiness of a given solution starting from a given machine end time table and a starting index */
/** */
long int PfspInstance::computeWT(std::vector< int > & sol, std::vector< std::vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i)
{    
   int j,m;
   long int wt;
   int jobNumber;

  // std::vector< std::vector < int >> previousMachineEndTimeMatrix(previousMachineEndTimeMatri);
   int prevj = previousMachineEndTimeMatrix[1][start_i-1];
   for(j=start_i;j<end_i;j++)
   {
       jobNumber = sol[j];
       prevj = prevj + processingTimesMatrix[jobNumber][1];
       previousMachineEndTimeMatrix[1][j] = prevj;
   }

     for ( j = start_i; j <= nbJob; ++j )
       {
           long int previousJobEndTime = previousMachineEndTimeMatrix[1][j];

           jobNumber = sol[j];
           for ( m = 2; m <= nbMac; ++m )
           {


           if ( previousMachineEndTimeMatrix[m][j-1] > previousJobEndTime )
           {
               previousMachineEndTimeMatrix[m][j] = previousMachineEndTimeMatrix[m][j-1] + processingTimesMatrix[jobNumber][m];

           }
           else
           {
               previousMachineEndTimeMatrix[m][j] = previousJobEndTime + processingTimesMatrix[jobNumber][m];
           }
           previousJobEndTime = previousMachineEndTimeMatrix[m][j];
       }
   }

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTimeMatrix[nbMac][j] - dueDates[sol[j]], 0L) * priority[sol[j]]);


    return wt;
}

void PfspInstance::computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail)
{
    int j,m;

    int jobNumber;
    int end_i = nbJob-1;
   // std::vector< std::vector < int >> head(previousMachineEndTimeMatri);
    int prevj = 0;
    int postj = 0;
    int k;
    for(j=1;j<nbJob;j++)
    {
        k = nbJob-j;
        jobNumber = sol[j];
        prevj = prevj + processingTimesMatrix[jobNumber][1];
        postj = postj + processingTimesMatrix[sol[k]][nbMac];
        head[1][j] = prevj;
        tail[nbMac][k] = postj;
    }

      for ( j = 1; j <= end_i; ++j )
        {
            k = nbJob-j;
            long int previousJobEndTime = head[1][j];
            long int postJobEndTime = tail[nbMac][k];

            jobNumber = sol[j];

            for ( m = 2; m <= nbMac; ++m )
            {
                int n = nbMac-m+1;

            if ( head[m][j-1] > previousJobEndTime )
            {
                head[m][j] = head[m][j-1] + processingTimesMatrix[jobNumber][m];

            }
            else
            {
                head[m][j] = previousJobEndTime + processingTimesMatrix[jobNumber][m];
            }

            if ( tail[n][k+1] > postJobEndTime )
            {
                tail[n][k] = tail[n][k+1] + processingTimesMatrix[sol[k]][n];

            }
            else
            {
                tail[n][k] = postJobEndTime + processingTimesMatrix[sol[k]][n];
            }

            previousJobEndTime = head[m][j];
            postJobEndTime = tail[n][k];
        }
    }
}

void PfspInstance::computeTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size)
{
    int j,m;

    int jobNumber;
    int end_i = size;
   // std::vector< std::vector < int >> head(previousMachineEndTimeMatri);
    int prevj = 0;
    int postj = 0;
    int k;
    for(j=1;j<size;j++)
    {
        k = size-j;
        jobNumber = sol[j];
        prevj = prevj + processingTimesMatrix[jobNumber][1];
        postj = postj + processingTimesMatrix[sol[k]][nbMac];
        head[1][j] = prevj;
        tail[nbMac][k] = postj;
    }

      for ( j = 1; j < end_i; ++j )
        {
            k = size-j;
            long int previousJobEndTime = head[1][j];
            long int postJobEndTime = tail[nbMac][k];

            jobNumber = sol[j];

            for ( m = 2; m <= nbMac; ++m )
            {
                int n = nbMac-m+1;
                if(k+1>=size)
                {
                    tail[n][k+1] = 0;
                }


            if ( head[m][j-1] > previousJobEndTime )
            {
                head[m][j] = head[m][j-1] + processingTimesMatrix[jobNumber][m];

            }
            else
            {
                head[m][j] = previousJobEndTime + processingTimesMatrix[jobNumber][m];
            }

            if ( tail[n][k+1] > postJobEndTime )
            {
                tail[n][k] = tail[n][k+1] + processingTimesMatrix[sol[k]][n];

            }
            else
            {
                tail[n][k] = postJobEndTime + processingTimesMatrix[sol[k]][n];
            }

            previousJobEndTime = head[m][j];
            postJobEndTime = tail[n][k];
        }
    }
}

void PfspInstance::computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail)
{
    computeNoIdleTAmatrices(sol,head,tail,nbJob);
}

void PfspInstance::computeNoIdleTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size)
{
    int j,m;

    int jobNumber;
    int end_i = size;//nbJob-1;
   // std::vector< std::vector < int >> head(previousMachineEndTimeMatri);
    long int a_h = 0;
    long int a_t = 0;
    int prevj = 0;
    int postj = 0;
    int k;
    for(j=1;j<size;j++)
    {
        k = size-j;
        jobNumber = sol[j];
        prevj = prevj + processingTimesMatrix[jobNumber][1];
        postj = postj + processingTimesMatrix[sol[k]][nbMac];
        head[1][j] = prevj;
        tail[nbMac][k] = postj;
    }

    k = size-1;//nbJob-1;



       for ( m = 2; m <= nbMac; ++m )
       {
           int n = nbMac-m+1;
           head[m][1] = head[m-1][1] + processingTimesMatrix[sol[1]][m];
           tail[n][k] = tail[n+1][k] + processingTimesMatrix[sol[k]][n];
       }

      for ( j = 2; j <= end_i; ++j )
        {
            k = size-j;
            long int previousJobEndTime = head[1][j];
            long int postJobEndTime = tail[nbMac][k];
            a_t = 0;
            a_h = 0;
            jobNumber = sol[j];

            for ( m = 2; m <= nbMac; ++m )
            {
                int n = nbMac-m+1;

            if ( head[m][j-1]+a_h > previousJobEndTime )
            {
                head[m][j] = head[m][j-1] + processingTimesMatrix[jobNumber][m]+ a_h;
            }
            else
            {
                head[m][j] = previousJobEndTime + processingTimesMatrix[jobNumber][m];
                a_h += previousJobEndTime-(head[m][j-1]+a_h);
            }

            if ( tail[n][k+1]+a_t > postJobEndTime )
            {
                tail[n][k] = tail[n][k+1] + processingTimesMatrix[sol[k]][n]+a_t;
            }
            else
            {
                tail[n][k] = postJobEndTime + processingTimesMatrix[sol[k]][n];
                a_t += postJobEndTime - (tail[n][k+1]+a_t);
            }

            previousJobEndTime = head[m][j];
            postJobEndTime = tail[n][k];
        }
    }

}

void PfspInstance:: computeSDSTTAmatrices(std::vector<int> &sol,std::vector< std::vector < int > >& head, std::vector< std::vector< int > >& tail,int size)
{
    int j,m;

    int jobNumber;
    int end_i = size;
   // std::vector< std::vector < int >> head(previousMachineEndTimeMatri);
    int prevj = 0;
    int postj = 0;
    int k,kp1;
    for(j=1;j<size;j++)
    {
        k = size-j;
        kp1 = k+1==size?0:k+1;
        jobNumber = sol[j];
        prevj = prevj + processingTimesMatrix[jobNumber][1]+ setUpTimes[1][sol[j-1]][jobNumber];
        postj = postj + processingTimesMatrix[sol[k]][nbMac] + setUpTimes[nbMac][sol[k]][sol[kp1]];
        head[1][j] = prevj;
        tail[nbMac][k] = postj;
    }

      for ( j = 1; j < end_i; ++j )
        {
            k = size-j;
            kp1 = k+1==size?0:k+1;
            long int previousJobEndTime = head[1][j];
            long int postJobEndTime = tail[nbMac][k];

            jobNumber = sol[j];

            for ( m = 2; m <= nbMac; ++m )
            {
                int n = nbMac-m+1;
                if(k+1>=size)
                {
                    tail[n][k+1] = 0;
                }

            long stpluspme = setUpTimes[m][sol[j-1]][jobNumber] + head[m][j-1];
            if ( previousJobEndTime > stpluspme )
            {
                head[m][j] = previousJobEndTime + processingTimesMatrix[jobNumber][m];

            }
            else
            {
                head[m][j] = stpluspme + processingTimesMatrix[jobNumber][m];
            }

            stpluspme = setUpTimes[n][sol[k]][sol[kp1]] + tail[n][kp1] ;
            if ( stpluspme > postJobEndTime )
            {
                tail[n][k] = stpluspme + processingTimesMatrix[sol[k]][n];

            }
            else
            {
                tail[n][k] = postJobEndTime + processingTimesMatrix[sol[k]][n];
            }

            previousJobEndTime = head[m][j];
            postJobEndTime = tail[n][k];
        }
    }
}

void PfspInstance::computeSDSThead(std::vector<int> &sol,std::vector< std::vector < int > >& head, int size)
{
    int j,m;
    int jobNumber;
    int end_i = size;
    int prevj = 0;
    for(j=1;j<size;j++)
    {
        jobNumber = sol[j];
        prevj = prevj + processingTimesMatrix[jobNumber][1]+ setUpTimes[1][sol[j-1]][jobNumber];
        head[1][j] = prevj;
    }

    for ( j = 1; j < end_i; ++j )
    {

        long int previousJobEndTime = head[1][j];
        jobNumber = sol[j];

        for ( m = 2; m <= nbMac; ++m )
        {
            long stpluspme = setUpTimes[m][sol[j-1]][jobNumber] + head[m][j-1];
            if ( previousJobEndTime > stpluspme )
            {
                head[m][j] = previousJobEndTime + processingTimesMatrix[jobNumber][m];

            }
            else
            {
                head[m][j] = stpluspme + processingTimesMatrix[jobNumber][m];
            }
            previousJobEndTime = head[m][j];
        }
    }
}


void inline computeTailss(std::vector<int> &sol, int size,std::vector< std::vector< int > > & tail,std::vector< std::vector< long> >& processingTimesMatrix,int nbMac)
{
    int j,m;

    int jobNumber;

    int postj = 0;

    for(j=size;j>=1;j--)
    {

        jobNumber = sol[j];
        postj = postj + processingTimesMatrix[jobNumber][nbMac];
        tail[nbMac][j] = postj;
    }

      for ( j = size; j >=1; --j )
        {


            long int postJobEndTime = tail[nbMac][j];

            jobNumber = sol[j];

            for ( m = nbMac-1; m >= 1; --m )
            {

            if ( tail[m][j+1] > postJobEndTime )
            {
                tail[m][j] = tail[m][j+1] + processingTimesMatrix[jobNumber][m];

            }
            else
            {
                tail[m][j] = postJobEndTime + processingTimesMatrix[jobNumber][m];
            }
            postJobEndTime = tail[m][j];
        }
    }
}

void PfspInstance::computeTails(std::vector<int> &sol, std::vector<std::vector<std::vector<int> > > &tails)
{
    for (int j = nbJob-1 ; j > 1; --j) {
        computeTailss(sol,j-1,tails[j],processingTimesMatrix,nbMac);
    }
}

/**
  // others machines :
  for ( m = 2; m <= nbMac; ++m )
  {
       previousJobEndTime = previousMachineEndTimeMatrix[m][start_i-1];

      for ( j = start_i; j <= nbJob; ++j )
      {
          jobNumber = sol[j];

          if ( previousMachineEndTimeMatrix[m-1][j] > previousJobEndTime )
          {
              previousMachineEndTimeMatrix[m][j] = previousMachineEndTimeMatrix[m-1][j] + processingTimesMatrix[jobNumber][m];
              previousJobEndTime = previousMachineEndTimeMatrix[m][j];
          }
          else
          {
              previousJobEndTime += processingTimesMatrix[jobNumber][m];
              previousMachineEndTimeMatrix[m][j] = previousJobEndTime;
          }
      }
  }*/


long int PfspInstance::computeWT(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime)
{
    int j, m;
    int jobNumber;
    long int wt;

    /**  And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /**  1st machine : */
    jobNumber = sol[job];
    previousMachineEndTime[job] = prevJob[1] + processingTimesMatrix[jobNumber][1];
    prevJob[1] = previousMachineEndTime[job];// -> qua iniziare ad aggiornare prevjob[machine 1]
    for ( j = job+1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];
    }


    /**  others machines : */
    for ( m = 2; m <= nbMac; ++m )
    {

        previousJobEndTime = prevJob[m]; // qua è uguale a prevjob[machine 2]
        jobNumber = sol[job];
        if ( previousMachineEndTime[job] > previousJobEndTime )
        {
            previousMachineEndTime[job] = previousMachineEndTime[job] + processingTimesMatrix[jobNumber][m];
            previousJobEndTime = previousMachineEndTime[job];
        }
        else
        {
            previousJobEndTime += processingTimesMatrix[jobNumber][m];
            previousMachineEndTime[job] = previousJobEndTime;
        }
        prevJob[m] = previousMachineEndTime[job];

        //j deve essere Job+1
        for ( j = job+1; j <= nbJob; ++j )
        {
            jobNumber = sol[j];

            if ( previousMachineEndTime[j] > previousJobEndTime )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime += processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }
        }
    }

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}

long int PfspInstance::computeWT(std::vector<int> &sol, std::vector<int>& makespans,int size)
{
    int j;
    long int wt=0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max( makespans[j] - dueDates[sol[j]] , 0L) * priority[sol[j]]);
    }

    return wt;
}

//this function sets up the prevJob and previousMachineEndTimestd::vector so that they can be used by the function above
void PfspInstance::computeWTs(std::vector<int> &sol,std::vector<int>& prevJob,int job,std::vector<int>& previousMachineEndTime)
{
    int j, m;
    int jobNumber;

    /**  And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /**  1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= job; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];                           
    }
     prevJob[1]= previousMachineEndTime[job];

    /**  others machines : */
    for ( m = 2; m <= nbMac; ++m )
    {
        previousMachineEndTime[1] +=
                processingTimesMatrix[sol[1]][m];
        previousJobEndTime = previousMachineEndTime[1];

        for ( j = 2; j <= job; ++j )
        {
            jobNumber = sol[j];

            if ( previousMachineEndTime[j] > previousJobEndTime )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime += processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }                                   
        }
            prevJob[m]= previousMachineEndTime[job];
    }

}

int PfspInstance::computeIdleTimeCoeff(std::vector<int>& prevJob, int job)
{
    int alpha = 0;

    int previousJobEndTime = prevJob[1] + processingTimesMatrix[job][1];

    for (int m = 2; m <= nbMac; ++m )
    {

            if ( prevJob[m] > previousJobEndTime )
            {

                previousJobEndTime = prevJob[m] + processingTimesMatrix[job][m];
                alpha++;
            }
            else
            {
                previousJobEndTime += processingTimesMatrix[job][m];

            }

    }
    return alpha;

}
/** Compute flowtime*/
long int PfspInstance::computeFT(std::vector< int >& sol)
{
    /** TO MODIFY IF WE HAVE TO DEAL WITH RELEASE DATES!!!*/
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += previousMachineEndTime[j];

    return wt;
}

long int PfspInstance::computeFT(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;
    for ( j = 1; j<= size; ++j )
        wt += previousMachineEndTime[j];

    return wt;
}

/**  Compute the weighted completion time of a given solution */
long int PfspInstance::computeWCT(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (previousMachineEndTime[j]  * priority[sol[j]]);

    return wt;
}

/** compute partial weighted completion time*/
long int PfspInstance::computeWCT(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (previousMachineEndTime[j]  * priority[sol[j]]);
    }

    return wt;

}

long int PfspInstance::computeWCT(std::vector<int> &sol, std::vector<int>& makespans,int size)
{
    int j;
    long int wt=0;

    for ( j = 1; j<= size; ++j ){

        wt +=  makespans[j] * priority[sol[j]];
    }

    return wt;
}

/**  total completion time*/
long int PfspInstance::computeTCT(std::vector< int > &sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;

    for ( j = 1; j<= nbJob; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

long int PfspInstance::computeTCT(std::vector< int > &sol,int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<=size; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

long int PfspInstance::computeTCT(std::vector<int> &sol, std::vector<int>& makespans,int size)
{
    int j;
    long int wt=0;

    for ( j = 1; j<=size; ++j ){

        wt += makespans[j];
    }

    return wt;
}


/**  Compute the weighted earliness of a given solution */
long int PfspInstance::computeWE(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);


    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);

    return wt;
}

/** compute partial weighted tardiness*/
long int PfspInstance::computeWE(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);


    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);
    }

    return wt;

}

long int PfspInstance::computeWE(std::vector<int> &sol, std::vector<int>& makespans,int size)
{
    int j;
    long int wt=0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - makespans[j] , 0L) * priority[sol[j]]);
    }

    return wt;
}

/**  Compute the weighted tardiness of a given solution */
long int PfspInstance::computeT(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L));// * priority[sol[j]]);

    return wt;
}

/** compute partial  tardiness*/
long int PfspInstance::computeT(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) );//**  priority[sol[j]]);
    }

    return wt;

}

long int PfspInstance::computeT(std::vector<int> &sol, std::vector<int>& makespans,int size)
{
    int j;
    long int wt=0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(makespans[j] - dueDates[sol[j]], 0L) );//**  priority[sol[j]]);
    }

    return wt;
}



/**  Compute the earliness of a given solution */
long int PfspInstance::computeE(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );

    return wt;
}

/** compute partial earliness */
long int PfspInstance::computeE(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );
    }

    return wt;

}


long int PfspInstance::computeE(std::vector<int> &sol, std::vector<int>& makespans,int size)
{
    int j;
    long int wt=0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]]-makespans[j], 0L) );//**  priority[sol[j]]);
    }

    return wt;
}

/** No wait permutation flowshop*/

/**  compute partials completition time distances*/

inline void computeNoWaitTimeDistances(std::vector<int> &sol, int nbMac, int size,std::vector<std::vector < long int > >& processingTimesMatrix,std::vector< int > & completionTimeDistance)
{
    for(int m=1; m <= nbMac;m++)
    {
        completionTimeDistance[1] += processingTimesMatrix[sol[1]][m];
    }

    for(int j=2; j <= size ; j++)
    {
        int max_ctd = 0;
        int i = sol[j-1];
        int jj = sol[j];
        for (int k = 1; k<= nbMac ; k++)
        {
            int temp_ctd = processingTimesMatrix[i][k];
            for(int h=k; h <= nbMac ; h++ )
            {
                temp_ctd += processingTimesMatrix[jj][h]-processingTimesMatrix[i][h];
            }
            if(temp_ctd > max_ctd){
                max_ctd = temp_ctd;
            }
        }
        completionTimeDistance[j] = max_ctd;
    }
}

/**  makespan */

long int PfspInstance::computeNWMS(std::vector<int> &sol)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,nbJob,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    for( int j=1; j<= nbJob ; j++)
    {
        nwms += completionTimeDistance[j];
    }

    return nwms;
}

/** partial makespan*/

long int PfspInstance::computeNWMS(std::vector<int> &sol, int size)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,size,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    for( int j=1; j<= size ; j++)
    {
        nwms += completionTimeDistance[j];
    }

    return nwms;
}

/**  No wait weighted tardiness*/
long int PfspInstance::computeNWWT(std::vector< int > &sol)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,nbJob,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=nbJob ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += (std::max(nwms - dueDates[sol[j]], 0L) * priority[sol[j]]);
    }
    return wt;
}

long int PfspInstance::computeNWWT(std::vector< int > &sol, int size)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,size,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=size ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += (std::max(nwms - dueDates[sol[j]], 0L) * priority[sol[j]]);
    }
    return wt;
}

/** No wait weighted earliness*/
long int PfspInstance::computeNWWE(std::vector< int > & sol)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,nbJob,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=nbJob ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += (std::max(dueDates[sol[j]] - nwms , 0L) * priority[sol[j]]);
    }
    return wt;
}

long int PfspInstance::computeNWWE(std::vector<int> &sol, int size)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,size,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=size ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += (std::max(dueDates[sol[j]] - nwms , 0L) * priority[sol[j]]);
    }

    return wt;
}

/**  No wait earliness*/

long int PfspInstance::computeNWE(std::vector< int > & sol)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,nbJob,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=nbJob ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += std::max(dueDates[sol[j]] - nwms , 0L);
    }

    return wt;
}

long int PfspInstance::computeNWE(std::vector<int> &sol, int size)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,size,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=size ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += std::max(dueDates[sol[j]] - nwms , 0L);
    }

    return wt;
}

/** No wait tardiness*/

long int PfspInstance::computeNWT(std::vector<int> &sol)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,nbJob,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=nbJob ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += (std::max(nwms - dueDates[sol[j]], 0L) );
    }

    return wt;
}

long int PfspInstance::computeNWT(std::vector<int> &sol, int size)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,size,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=size ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += (std::max(nwms - dueDates[sol[j]], 0L));
    }

    return wt;
}

/**  total completion time*/
long int PfspInstance::computeNWTCT(std::vector< int > &sol)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,nbJob,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=nbJob ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += nwms;
    }

    return wt;
}

long int PfspInstance::computeNWTCT(std::vector< int > &sol,int size)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,size,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=size ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += nwms;
    }

    return wt;
}

long int PfspInstance::computeNWWCT(std::vector< int > &sol)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,nbJob,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=nbJob ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += nwms * priority[sol[j]];
    }

    return wt;
}

long int PfspInstance::computeNWWCT(std::vector< int > &sol,int size)
{
    std::vector < int > completionTimeDistance(nbJob+1,0);

    computeNoWaitTimeDistances(sol,nbMac,size,processingTimesMatrix,completionTimeDistance);

    int nwms = 0;
    int wt = 0;

    for(int j=1 ; j<=size ; j++ )
    {
        nwms += completionTimeDistance[j];
        wt += nwms* priority[sol[j]];
    }

    return wt;
}

/** No idle permutation flowshop*/


inline std::vector< long int > computeNoIdlePartialMakespans(std::vector< int >& sol,std::vector< std::vector <long int> >& processingTimesMatrix,int nbMac,int nbJob,int size)
{
    std::vector< long int > partialMs(nbJob+1,0);

    long int nims = processingTimesMatrix[sol[1]][1];
    partialMs[1] = nims;
    std::vector< int > minimumDiff(nbMac,0);

    for(int m=1;m<nbMac;++m)
    {
        minimumDiff[m] = processingTimesMatrix[sol[1]][m+1];
        partialMs[1] += processingTimesMatrix[sol[1]][m+1];
    }


    for(int j=2; j<= size; ++j)
    {
        nims += processingTimesMatrix[sol[j]][1];
        partialMs[j] = nims;

        for(int m=1;m<nbMac;++m)
        {
            minimumDiff[m] = std::max(minimumDiff[m]-processingTimesMatrix[sol[j]][m],0L) + processingTimesMatrix[sol[j]][m+1];
            partialMs[j] += minimumDiff[m];
        }
    }

    /** for(int m=1;m<nbMac;++m)
        nims += minimumDiff[m];

    return nims;*/
    return partialMs;
}


inline std::vector< long int > computeNoIdlePartialMakespans(std::vector< int >& sol,std::vector< std::vector <long int> >& processingTimesMatrix,int nbMac,int nbJob)
{
    return computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob,nbJob);
}


long int PfspInstance::computeNIMS(std::vector<int> &sol)
{
   long int nims = processingTimesMatrix[sol[1]][1];

   std::vector< int > minimumDiff(nbMac,0);
   for(int m=1;m<nbMac;++m)
       minimumDiff[m] = processingTimesMatrix[sol[1]][m+1];


   for(int j=2; j<= nbJob; ++j)
   {
       nims += processingTimesMatrix[sol[j]][1];
       for(int m=1;m<nbMac;++m)
       {
           minimumDiff[m] = std::max(minimumDiff[m]-processingTimesMatrix[sol[j]][m],0L) + processingTimesMatrix[sol[j]][m+1];

       }
   }

   for(int m=1;m<nbMac;++m)
       nims += minimumDiff[m];

   return nims;
}

long int PfspInstance::computeNIMS(std::vector<int> &sol,int size)
{
   long int nims = processingTimesMatrix[sol[1]][1];

   std::vector< int > minimumDiff(nbMac,0);
   for(int m=1;m<nbMac;++m)
       minimumDiff[m] = processingTimesMatrix[sol[1]][m+1];


   for(int j=2; j<= size; ++j)
   {
       nims += processingTimesMatrix[sol[j]][1];
       for(int m=1;m<nbMac;++m)
       {
           minimumDiff[m] = std::max(minimumDiff[m]-processingTimesMatrix[sol[j]][m],0L) + processingTimesMatrix[sol[j]][m+1];
       }
   }

   for(int m=1;m<nbMac;++m)
       nims += minimumDiff[m];

   return nims;
}



long int PfspInstance::computeNIMS(std::vector<int> &sol, long int nims)
{
    std::vector< int > minimumDiff(nbMac,0);
    for(int m=1;m<nbMac;++m)
        minimumDiff[m] = processingTimesMatrix[sol[1]][m+1];


    for(int j=2; j<= nbJob; ++j)
    {
        for(int m=1;m<nbMac;++m)
        {
            minimumDiff[m] = std::max(minimumDiff[m]-processingTimesMatrix[sol[j]][m],0L) + processingTimesMatrix[sol[j]][m+1];

        }
    }

    for(int m=1;m<nbMac;++m)
        nims += minimumDiff[m];    


    return nims;
}


// Compute no idle weighted tardiness
long int PfspInstance::computeNIWT(std::vector<int> &sol)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob);
    long int wt = 0;
    for ( int j = 1; j<= nbJob; ++j )
        wt += (std::max(partials[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}

long int PfspInstance::computeNIWT(std::vector<int> &sol, int size)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob,size);
    long int wt = 0;
    for ( int j = 1; j<= size; ++j )
        wt += (std::max(partials[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}

// Compute no idle weighted earliness
long int PfspInstance::computeNIWE(std::vector< int > & sol)
{
    std::vector< long int > previousMachineEndTime = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob);
    long int wt = 0;
    for (int j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);

    return wt;
}

long int PfspInstance::computeNIWE(std::vector<int> &sol, int size)
{
    std::vector< long int > previousMachineEndTime = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob,size);
    long int wt = 0;
    for (int j = 1; j<= size; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);

    return wt;
}

// Compute no idle earliness

long int PfspInstance::computeNIE(std::vector< int > & sol)
{
    std::vector< long int > previousMachineEndTime = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob);
    long int wt = 0;
    for (int j = 1; j<= nbJob; ++j )
        wt += std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L);

    return wt;
}

long int PfspInstance::computeNIE(std::vector<int> &sol, int size)
{
    std::vector< long int > previousMachineEndTime = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob,size);
    long int wt = 0;
    for (int j = 1; j<= size; ++j )
        wt += std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L);

    return wt;
}

// Compute no idle tardiness

long int PfspInstance::computeNIT(std::vector<int> &sol)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob);
    long int wt = 0;
    for ( int j = 1; j<= nbJob; ++j )
        wt += (std::max(partials[j] - dueDates[sol[j]], 0L) );

    return wt;
}

long int PfspInstance::computeNIT(std::vector<int> &sol, int size)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob,size);
    long int wt = 0;
    for ( int j = 1; j<= size; ++j )
        wt += (std::max(partials[j] - dueDates[sol[j]], 0L));

    return wt;
}

// Compute no idle total completion time

long int PfspInstance::computeNITCT(std::vector<int> &sol)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob);
    long int wt = 0;
    for ( int j = 1; j<= nbJob; ++j )
        wt += partials[j] ;

    return wt;
}

long int PfspInstance::computeNITCT(std::vector<int> &sol, int size)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob,size);
    long int wt = 0;
    for ( int j = 1; j<= size; ++j )
        wt += partials[j];

    return wt;
}

long int PfspInstance::computeNIWCT(std::vector<int> &sol)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob);
    long int wt = 0;
    for ( int j = 1; j<= nbJob; ++j )
        wt += (partials[j] * priority[sol[j]]) ;

    return wt;
}

long int PfspInstance::computeNIWCT(std::vector<int> &sol, int size)
{
    std::vector< long int > partials = computeNoIdlePartialMakespans(sol,processingTimesMatrix,nbMac,nbJob,size);
    long int wt = 0;
    for ( int j = 1; j<= size; ++j )
        wt += (partials[j] * priority[sol[j]]);

    return wt;
}

// Compute sequence dependent setup times

long int PfspInstance::computeSDSTMS(std::vector<int> &sol)
{
    int j, m;
    int jobNumber;

    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /**  1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1] + setUpTimes[1][sol[j-1]][jobNumber];
    }

    /**  others machines : */
    for ( m = 2; m <= nbMac; ++m )
    {
        previousMachineEndTime[1] += processingTimesMatrix[sol[1]][m];
        previousJobEndTime = previousMachineEndTime[1];


        for ( j = 2; j <= nbJob; ++j )
        {
            jobNumber = sol[j];
            long stpluspme = setUpTimes[m][sol[j-1]][jobNumber] + previousJobEndTime;
            if (previousMachineEndTime[j] > stpluspme )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime = stpluspme + processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }
        }
    }

    return previousMachineEndTime[nbJob];
}

long int PfspInstance::computeSDSTMS(std::vector<int> &sol, int size)
{
    int j, m;
    int jobNumber;

    // We need end times on previous machine :
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    // And the end time of the previous job, on the same machine :
    long int previousJobEndTime;
     //1st machine :

    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= size; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1] + setUpTimes[1][sol[j-1]][jobNumber];
    }

    // others machines :
    for ( m = 2; m <= nbMac; ++m )
    {
        previousMachineEndTime[1] += processingTimesMatrix[sol[1]][m];
        previousJobEndTime = previousMachineEndTime[1];


        for ( j = 2; j <= size; ++j )
        {
            jobNumber = sol[j];

            long stpluspme = setUpTimes[m][sol[j-1]][jobNumber] + previousJobEndTime;
            if (previousMachineEndTime[j] > stpluspme )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime = stpluspme + processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }
        }
    }
    return previousMachineEndTime[size];
}

inline void computePartialSDSTMakespans(std::vector< int >& sol,std::vector< long int >& previousMachineEndTime,std::vector<std::vector< long> >& processingTimesMatrix,std::vector< std::vector< std::vector < int > > >& setUpTimes,int nbJob, int nbMac)
{
    int j, m;
    int jobNumber;
   /**  And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /**  1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1] + setUpTimes[1][sol[j-1]][jobNumber];
    }

    /**  others machines : */
    for ( m = 2; m <= nbMac; ++m )
    {
        previousMachineEndTime[1] += processingTimesMatrix[sol[1]][m];
        previousJobEndTime = previousMachineEndTime[1];


        for ( j = 2; j <= nbJob; ++j )
        {
            jobNumber = sol[j];
            long stpluspme = setUpTimes[m][sol[j-1]][jobNumber] + previousJobEndTime;
            if (previousMachineEndTime[j] > stpluspme )
            {
                previousMachineEndTime[j] = previousMachineEndTime[j] + processingTimesMatrix[jobNumber][m];
                previousJobEndTime = previousMachineEndTime[j];
            }
            else
            {
                previousJobEndTime = stpluspme + processingTimesMatrix[jobNumber][m];
                previousMachineEndTime[j] = previousJobEndTime;
            }
        }
    }
}

/**  Compute the weighted completion time of a given solution */
long int PfspInstance::computeSDSTWCT(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (previousMachineEndTime[j]  * priority[sol[j]]);

    return wt;
}

/** compute partial weighted completion time*/
long int PfspInstance::computeSDSTWCT(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (previousMachineEndTime[j]  * priority[sol[j]]);
    }

    return wt;

}

/**  total completion time*/
long int PfspInstance::computeSDSTTCT(std::vector< int > &sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;

    for ( j = 1; j<= nbJob; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

long int PfspInstance::computeSDSTTCT(std::vector< int > &sol,int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<=size; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

/**  Compute the weighted earliness of a given solution */
long int PfspInstance::computeSDSTWE(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);


    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);

    return wt;
}

/** compute partial weighted tardiness*/
long int PfspInstance::computeSDSTWE(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);


    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);
    }

    return wt;

}

/**  Compute the weighted tardiness of a given solution */
long int PfspInstance::computeSDSTT(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L));// * priority[sol[j]]);

    return wt;
}

/** compute partial  tardiness*/
long int PfspInstance::computeSDSTT(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) );//**  priority[sol[j]]);
    }

    return wt;

}

/**  Compute the earliness of a given solution */
long int PfspInstance::computeSDSTE(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );

    return wt;
}

/** compute partial earliness */
long int PfspInstance::computeSDSTE(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );
    }

    return wt;

}

/**  Compute the weighted tardiness of a given solution */
long int PfspInstance::computeSDSTWT(std::vector< int > & sol)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 );
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}

/** compute partial weighted tardiness*/
long int PfspInstance::computeSDSTWT(std::vector<int> &sol, int size)
{
    int j;
    long int wt;
    /**  We need end times on previous machine : */
   std::vector< long int > previousMachineEndTime ( nbJob + 1 ,0);
    /**  And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);


    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);
    }

    return wt;

}

long int PfspInstance::computeMSLB()
{
    int max_j = 0;
    int max_i[nbMac];

    for(int i=0;i<nbJob;i++)
    {
        int temp_j=0;
        for(int j=1;j<=nbMac;j++)
        {
            int pt = processingTimesMatrix[i][j];
            max_i[j-1]=i==0?pt:pt+max_i[j-1];
            temp_j+=pt;
        }
        max_j = temp_j>max_j?temp_j:max_j;
    }
    int maxx_i = 0;
    for(int i=0;i<nbMac;i++)
    {
        maxx_i = max_i[i]>maxx_i?max_i[i]:maxx_i;
    }
    int P=max_j>maxx_i?max_j:maxx_i;
    return P;
}

void PfspInstance::computeHead(std::vector<int>& sol,std::vector< std::vector< int > >& head, int njobs)
{
    int j,m;
    int jobNumber;
    int prevj = 0;
    for(j=1;j<njobs;j++)
    {
        jobNumber = sol[j];
        prevj = prevj + processingTimesMatrix[jobNumber][1];
        head[1][j] = prevj;
    }

      for ( j = 1; j < njobs; ++j )
        {
            long int previousJobEndTime = head[1][j];
            jobNumber = sol[j];

            for ( m = 2; m <= nbMac; ++m )
            {
                if ( head[m][j-1] > previousJobEndTime )
                {
                    head[m][j] = head[m][j-1] + processingTimesMatrix[jobNumber][m];
                }
                else
                {
                    head[m][j] = previousJobEndTime + processingTimesMatrix[jobNumber][m];
                }
                previousJobEndTime = head[m][j];
            }
    }
}




