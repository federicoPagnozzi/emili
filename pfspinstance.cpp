//#define RELDUE_FIELD
//#define TWO_VECTORS
//#define HORIZONTAL
#define IDLE_ONE_STEP
#define IDLE_TWO_STEPS
#define IDLE_FORWARD
#define IDLE_BACKWARDS
#define SIMPLE_IDLE_INSERTION

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


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include "pfspinstance.h"
#include <numeric>

// #define ENABLE_SSE 10

#ifdef ENABLE_SSE
#ifdef __SSE__
#include <xmmintrin.h>
#include <emmintrin.h>
#endif
#endif

using namespace std;

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

#pragma region NEWCODE

int PfspInstance::getNbStages()
{
	return nbStages;
}

void PfspInstance::setNbStages(int stagesCount)
{
	this->nbStages = stagesCount;
}

vector<int>& PfspInstance::getStages()
{
	return stages;
}

void PfspInstance::setStages(vector< int > &stages)
{
	this->stages = stages;
}

std::string newToLower(string taget)
{
	string temp = taget;
	std::transform(temp.begin(), temp.end(), temp.begin(), ::tolower);
	return temp;
}

std::string printVector(vector< int >& vec)
{
	std::string result = "";
	//for (int i = 1; i < vec.size(); i++) result += std::to_string(vec[i]) + " ";

	for (std::vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
		result += std::to_string(*it) + " ";

	return result;
}

std::string printVector(vector< long >& vec)
{
	std::string result = "";
	//for (int i = 1; i < vec.size(); i++) result += std::to_string(vec[i]) + " ";
	
	for (std::vector<long>::iterator it = vec.begin(); it != vec.end(); ++it)
		result += std::to_string(*it) + " ";

	return result;
}

std::string printJaggedVector(vector< vector < int > >& jaggedVector)
{
	string temp = "";

	for (std::vector<vector<int>>::iterator row = jaggedVector.begin(); row != jaggedVector.end(); ++row) {
		for (std::vector<int>::iterator col = row->begin(); col != row->end(); ++col) 
			temp += std::to_string(*col) + " ";
		temp += "\n";
	}

	/*for (std::vector<vector<int>>::iterator it = jaggedVector.begin(); it != jaggedVector.end(); ++it)
		for (std::vector<int>::iterator it = it.begin(); it != jaggedVector.end(); ++it)
			result += std::to_string(*it) + " ";
			*/

	/*for (int x = 0; x < dim1; x++)
	{
		temp += std::to_string(x) + "-> ";
		for (int y = 0; y < dimensions[x]; y++)
		{
			temp += std::to_string(jaggedVector[x][y]) + ", ";
		}
		temp += "\n";
	}*/
	return temp;
}

std::string printPairedVector(vector< pair< int, int > >& vec)
{
	string temp = "";

	for (std::vector< pair< int, int > >::iterator it = vec.begin(); it != vec.end(); ++it)
		temp += std::to_string(it-vec.begin()) + "(" + std::to_string(it->first) + "," + std::to_string(it->second) + "), ";

	/*for (int x = 0; x < dim1; x++)
		temp += "J" + std::to_string(x) + " " + std::to_string(vec[x].first) + "FT: " + std::to_string(vec[x].first) + " ;; ";*/
	return temp;
}
#pragma endregion NEWCODE

void PfspInstance::setSilence(bool s)
{
    this->silence = s;
}

/* Allow the memory for the processing times matrix : */
void PfspInstance::allowMatrixMemory(int nbJ, int nbM)
{
  processingTimesMatrix.resize(nbJ+1);          // 1st dimension

  for (int cpt = 0; cpt < nbJ+1; ++cpt)
    processingTimesMatrix[cpt].resize(nbM+1); // 2nd dimension

	dueDates.resize(nbJ+1);
	priority.resize(nbJ+1);
#pragma region NEWCODE
	releaseDates.resize(nbJ + 1);
	weightsE.resize(nbJ + 1);
	earlinessDD.resize(nbJ + 1);
	tardinessDD.resize(nbJ + 1);
#pragma endregion NEWCODE
}


long int PfspInstance::getTime(int job, int machine)
{
  if (job == 0)
    return 0;
  else
  {
    if ((job < 1) || (job > nbJob) || (machine < 1) || (machine > nbMac))
      cout    << "ERROR. file:pfspInstance.cpp, method:getTime. Out of bound. job=" << job
          << ", machine=" << machine << std::endl;

    return processingTimesMatrix[job][machine];
  }
}


/* Read the instance from file : */
bool PfspInstance::readDataFromFile(char * fileName)
{
	string code = "No Code";
	bool everythingOK = true;
	int j, m; // iterators
	long int readValue;
	string str;
	ifstream fileIn;

	char * aux2;
	char fileNameOK[100] = "";

	aux2 = (strrchr(fileName, '/'));

	if (aux2 == NULL)
		aux2 = fileName;
	else
		aux2 += 1;

	//strcat(fileNameOK, aux2);
#if defined(_WIN32) || defined(_WIN64)
	strcat_s(fileNameOK, aux2);
#else
        strcat(fileNameOK, aux2);
#endif
    if(!silence)
    {
	cout << "name : " << fileNameOK << endl;
	cout << "file : " << fileName << endl;
    }
	fileIn.open(fileName);

	if ( fileIn.is_open() ) {

#pragma region NEWCODE
		fileIn >> str;
		if (str == "HFSDDW" || str == "HFS")
		{
			code = str;
			fileIn >> nbJob;
		}
		else nbJob = stol(str);
		//fileIn >> nbJob;
#pragma endregion NEWCODE
        fileIn >> nbMac;
        fileIn >> readValue;
        if(readValue == 12345)
        {
            std::string fname(fileName);
            return readDataFromFile(fname);
        }        

#pragma region NEWCODE
		// We have stages  Jobs Machines Stages
		if (readValue != 0)
		{
			nbStages = readValue;

			vector< int > st(nbStages + 1);
			stages = st;
			
			fileIn >> readValue;// Just to read first index of the table.  OR next line of Stages.

			if (readValue != 0)
			{
				// We read the second line.
				stages[1] = readValue; 
				
				int sum = stages[1];
				for (int i = 2; i <= nbStages; i++)
				{
					fileIn >> stages[i];
					sum += stages[i];
				}
				//If machines in stages doesnt sum the amount of machines in the instance, instance is wrong. 
				if (sum != nbMac) throw new runtime_error("Instance machines and machines per stage doesn't fit");
				fileIn >> readValue;// Just to read first index of the table.  
			}
			else
			{
				// We have first line like 50  15   3   That means a HFS((PM^5)(PM^5)(PM^5)) 50 jobs, 15 machines distributed in 3 stages
				int machinesPerStage = nbMac / nbStages;

				// If machines is not divisible by stages, instance is wrong. 
				int mod = nbMac % nbStages;
				if (mod != 0) throw new runtime_error("Instance machines and machines per stage doesn't fit");

				for (int i = 1; i <= nbStages; i++)
					stages[i] = machinesPerStage;
			}
			stages[0] = 1; // Yes, a ghost stage to keep with Federicos practice. With only one machine. 
		}

		int nbRows = nbMac;
		if (nbStages > 0) nbRows = nbStages;
#pragma endregion NEWCODE
		allowMatrixMemory(nbJob, nbRows);
        if(!silence){
            cout << "File " << fileName << " is now open, start to read..." << std::endl;
            cout << "Number of jobs : " << nbJob << std::endl;
            cout << "Number of machines : " << nbMac << std::endl;
			cout << "Number of stages : " << nbStages << ". { " << printVector(stages) << " }" << std::endl;
            cout << "Memory allowed." << std::endl;
            cout << "Start to read matrix..." << std::endl;
        }
		for (j = 1; j <= nbJob; ++j)
		{
            for (m = 1; m <= nbRows; ++m)
			{
                if(!(j==1 && m==1))
                {
				fileIn >> readValue; // The number of each machine, not important !
                }
                fileIn >> readValue; // Process Time

				processingTimesMatrix[j][m] = readValue;
			}
		}
        fileIn >> str; 

#pragma region NEWCODE
		int LBCmax = -1;
		if ((newToLower(str).compare("lbcmax:") == 0))
		{
			fileIn >> LBCmax;
			fileIn >> str;
		}
		// While str is not Reldue and not end of file... Ignore everything... 
		while ((newToLower(str).compare("reldue") != 0) && (!fileIn.eof()))
		{
			fileIn >> str;
		}
		if (newToLower(str).compare("reldue") == 0)
		{
			for (j = 1; j <= nbJob; ++j)
			{
				fileIn >> readValue; // Release Date
				releaseDates[j] = readValue;
				fileIn >> readValue; // Due Date
				dueDates[j] = readValue;
				fileIn >> readValue; // Earliness Weights
				weightsE[j] = readValue;
				fileIn >> readValue; // Tardiness Weights
				priority[j] = readValue;
			}
		}

		// While str is not ddw and not end of file... Ignore everything... 
		while ((newToLower(str).compare("ddw") != 0) && (!fileIn.eof()))
		{
			fileIn >> str;
		}
		if (newToLower(str).compare("ddw") == 0)
		{
			for (j = 1; j <= nbJob; ++j)
			{
				fileIn >> readValue; // Earliness Due Date
				earlinessDD[j] = readValue;
				fileIn >> readValue; // Tardiness Due Date
				tardinessDD[j] = readValue;
			}
		}
#pragma endregion NEWCODE
        if(!silence)
        cout << "All is read from file." << std::endl;
		fileIn.close();
	}
	else
	{
        if(!silence)
		cout    << "ERROR. file:pfspInstance.cpp, method:readDataFromFile, "
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
    string str;
    ifstream fileIn;

    char * aux2;
    char fileNameOK[100] = "";

    aux2 = (strrchr(fileName, '/'));

    if (aux2 == NULL)
        aux2 = fileName;
    else
        aux2 += 1;

    //strcat_s(fileNameOK, aux2);
    //strcat(fileNameOK, aux2);
#if defined(_WIN32) || defined(_WIN64)
	strcat_s(fileNameOK, aux2);
#else
        strcat(fileNameOK, aux2);
#endif
    if(!silence)
    {
		cout << "name : " << fileNameOK << endl;
		cout << "file : " << fileName << endl;
    }
    fileIn.open(fileName);

    if ( fileIn.is_open() ) {

#pragma region NEWCODE
		fileIn >> str;
		if (str == "HFSDDW") fileIn >> nbJob;
		else nbJob = stol(str);
		//fileIn >> nbJob;
#pragma endregion NEWCODE
		fileIn >> nbMac;
		fileIn >> readValue;
		if (readValue == 12345)
		{
			std::string fname(fileName);
			return readDataFromFile(fname);
		}
#pragma region NEWCODE
		// We have stages  Jobs Machines Stages
		if (readValue != 0)
		{
			nbStages = readValue;
			vector< int > st(nbStages + 1);
			stages = st;

			fileIn >> readValue;// Just to read first index of the table.  OR next line of Stages.

			if (readValue != 0)
			{
				// We read the second line.
				stages[1] = readValue;

				int sum = stages[1];
				for (int i = 2; i <= nbStages; i++)
				{
					fileIn >> stages[i];
					sum += stages[i];
				}
				//If machines in stages doesnt sum the amount of machines in the instance, instance is wrong. 
				if (sum != nbMac) throw new runtime_error("Instance machines and machines per stage doesn't fit");
				fileIn >> readValue;// Just to read first index of the table.  
			}
			else
			{
				// We have first line like 50  15   3   That means a HFS((PM^5)(PM^5)(PM^5)) 50 jobs, 15 machines distributed in 3 stages
				int machinesPerStage = nbMac / nbStages;

				// If machines is not divisible by stages, instance is wrong. 
				int mod = nbMac % nbStages;
				if (mod != 0) throw new runtime_error("Instance machines and machines per stage doesn't fit");

				for (int i = 1; i <= nbStages; i++)
					stages[i] = machinesPerStage;
			}
			stages[0] = 1; // Yes, a ghost stage to keep with Federicos practice. With only one machine. 
		}

		int nbRows = nbMac;
		if (nbStages > 0) nbRows = nbStages;
#pragma endregion NEWCODE

        allowMatrixMemory(nbJob, nbRows);

        if(!silence){
            cout << "File " << fileName << " is now open, start to read..." << std::endl;
            cout << "Number of jobs : " << nbJob << std::endl;
			cout << "Number of machines : " << nbMac << std::endl;
			cout << "Number of stages : " << nbStages << ". { " << printVector(stages) << " }" << std::endl;
            cout << "Memory allowed." << std::endl;
            cout << "Start to read matrix..." << std::endl;
        }
        for (j = 1; j <= nbJob; ++j)
        {
			for (m = 1; m <= nbRows; ++m)
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
			setUpTimes.resize(nbRows + 1);
			for (int i = 0; i<nbRows + 1; i++){
                setUpTimes[i].resize(nbJob+1);
                for(int j=0;j<nbJob+1;j++)
                {
                    setUpTimes[i][j].resize(nbJob+1);
                }
            }

			for (int i = 1; i< nbRows + 1; i++)
            {
                fileIn >> str; // this is not read
                for(int j=1;j<nbJob+1;j++)
                {
                    for(int k=1;k<nbJob+1;k++)
                    {
                        fileIn >> readValue; // -1
                        //fileIn >> readValue;
                        setUpTimes[i][j][k] = readValue;
                    }
                }
            }

        }

		fileIn >> str;

#pragma region NEWCODE
		int LBCmax = -1;
		if ((newToLower(str).compare("lbcmax:") == 0))
		{
			fileIn >> LBCmax;
			fileIn >> str;
		}
		// While str is not Reldue and not end of file... Ignore everything... 
		while ((newToLower(str).compare("reldue") != 0) && (!fileIn.eof()))
		{
			fileIn >> str;
		}
		if (newToLower(str).compare("reldue") == 0)
		{
			for (j = 1; j <= nbJob; ++j)
			{
				fileIn >> readValue; // Release Date
				releaseDates[j] = readValue;
				fileIn >> readValue; // Due Date
				dueDates[j] = readValue;
				fileIn >> readValue; // Earliness Weights
				weightsE[j] = readValue;
				fileIn >> readValue; // Tardiness Weights
				priority[j] = readValue;
			}
		}

		// While str is not ddw and not end of file... Ignore everything... 
		while ((newToLower(str).compare("ddw") != 0) && (!fileIn.eof()))
		{
			fileIn >> str;
		}
		if (newToLower(str).compare("ddw") == 0)
		{
			for (j = 1; j <= nbJob; ++j)
			{
				fileIn >> readValue; // Earliness Due Date
				earlinessDD[j] = readValue;
				fileIn >> readValue; // Tardiness Due Date
				tardinessDD[j] = readValue;
			}
		}
#pragma endregion NEWCODE
        if(!silence)
        cout << "All is read from file." << std::endl;
        fileIn.close();
    }
    else
    {
        if(!silence)
        cout    << "ERROR. file:pfspInstance.cpp, method:readDataFromFile, "
                << "error while opening file " << fileName << std::endl;
        everythingOK = false;

    }

    return everythingOK;
}

// Weird exception to read some weird instances.
bool PfspInstance::readDataFromFile(const string _fileName)
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
          /*for (unsigned int j=1 ; j<nbMac ; j++)
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
    return true;
}


#ifndef ENABLE_SSE
// Compute Partial makespans returns previousMachineEndTime with the completion times of last machine of the flowshop.
inline void computePartialMakespans( vector< int >& sol, vector< long int >& previousMachineEndTime,vector< vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    long int previousJobEndTime;
    int j, m;
    int jobNumber;
    /* 1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];
    }

    /* others machines : */
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

//NEWCODE GOOD
/* Compute the weighted tardiness of a given solution */
long int PfspInstance::computeWT(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    /*!!! It could be implemented using SIMD instructions... */
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}
#else
inline void computePartialMakespans( vector< int >& sol, vector< long int >& previousMachineEndTime,vector< vector< long> >& processingTimesMatrix,int nbJob, int nbMac)
{
    long int previousJobEndTime;
    int j, m;
    int jobNumber;
    /* 1st machine : */
    previousMachineEndTime[0] = 0;
    int prevmj = 0;
    for ( j = 1; j <= nbJob; j++ )
    {
        prevmj += processingTimesMatrix[sol[j]][1];
        previousMachineEndTime[j] = prevmj;
    }

    /* others machines : */
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

/* Compute the weighted tardiness of a given solution */
long int PfspInstance::computeWT(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);
     wt = 0;
    /* USING FLOATS */
   __m128 a,b,z,p;
    z = _mm_setzero_ps();
    float res[4] __attribute__((aligned(16)));

    for (j=1; j<= nbJob ; j+=4)
    {
        //store 4 values in a,b,p
        a = _mm_set_ps(previousMachineEndTime[j],previousMachineEndTime[j+1],previousMachineEndTime[j+2],previousMachineEndTime[j+3]);
        b = _mm_set_ps(dueDates[sol[j]],dueDates[sol[j+1]],dueDates[sol[j+2]],dueDates[sol[j+3]]);
        p = _mm_set_ps(priority[sol[j]],priority[sol[j+1]],priority[sol[j+2]],priority[sol[j+3]]);
        // completion time - due date
        a = _mm_sub_ps(a,b);
        // max( tardiness , zero )
        a = _mm_max_ps(a,z);
        // tardiness * priority
        a = _mm_mul_ps(a,p);
        // add all 4 values
        a = _mm_add_ps(a, _mm_movehl_ps(a, a));
        a = _mm_add_ss(a, _mm_shuffle_ps(a, a, 1));
        // store on res
        _mm_store_ss(res,a);
        //final add
        wt += res[0];
    }

   /*
    __m128i a,b,p;
    __m128i z = _mm_setzero_si128();
    long int res[2] __attribute__((aligned(16)));
    for(j=1;j<=nbJob;j+=2)
    {
        a = _mm_loadu_si128(reinterpret_cast< __m128i*> (previousMachineEndTime.data()+j));
            b = _mm_loadu_si128(reinterpret_cast< __m128i*> (dueDates.data()+j));
            p = _mm_loadu_si128(reinterpret_cast< __m128i*> (priority.data()+j));

            a = _mm_sub_epi64(a, b);
            a = _mm_max_epi16(a, z);//This is not the right one and that's why it does not work!
            a = _mm_mul_epu32(a, p); // this maybe could do it...

            _mm_store_si128(reinterpret_cast<__m128i*> (res), a);
            wt+=res[0]+res[1];
    } */
    for (; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}
#endif

/*compute partial weighted tardiness*/
long int PfspInstance::computeWT(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 ,0);
    /* And the end time of the previous job, on the same machine : */
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

#pragma region NEWCODE

// Given an array obtains the index of min value 
int minValue(vector<int> &finishingTimes)
{
	//// This gives the min index 
	//int min_index = std::min_element(finishingTimes.begin() + 1, finishingTimes.end()) - finishingTimes.begin();

	int minValue = std::numeric_limits<int>::max();;
	int ID = 0;
	int size = finishingTimes.size();
	for (size_t i = 0; i < size; i++)
	{
		if (finishingTimes[i] < minValue)
		{
			minValue = finishingTimes[i];
			ID = i;
		}
	}
	return ID;
}

// Given an array obtains the index of max value 
int maxValue(vector<int> &finishingTimes)
{
	//// This gives the min index?
	//int max_index = std::max_element(finishingTimes.begin() + 1, finishingTimes.end()) - finishingTimes.begin();

	int maxValue = std::numeric_limits<int>::min();;
	int ID = 0;
	int size = finishingTimes.size();
	for (size_t i = 0; i < size; i++)
	{
		if (finishingTimes[i] > maxValue)
		{
			maxValue = finishingTimes[i];
			ID = i;
		}
	}
	return maxValue;
}

// Given an array and a valueCealing. Obtains the higher value that is under or equal to cealing, if there is none, will just return minimum value.
int maxLower(vector<int> &finishingTimes, int cealing)
{
	int minValue = std::numeric_limits<int>::max();
	int bestSuited = -1;
	int bestSuitedID = -1;
	int minID = 0;
	int size = finishingTimes.size();

	// Browse vector trying to get the ID of the maximum value under (or equal) the cealing  vec[2,4,6,8,10] cealing(7) -> 2 (id for 6 biggest element under 7)
	for (size_t i = 0; i < size; i++)
	{
		int value = finishingTimes[i];
		if (value < minValue)
		{
			minValue = value;
			minID = i;
		}
		if ((value <= cealing) && (value > bestSuited))
		{
			bestSuited = value;
			bestSuitedID = i;			
		}
	}

	if (bestSuitedID == -1) bestSuitedID = minID;

	return bestSuitedID;
}

// Given two arrays, orders both in ascending order based on the values of the second one.
void orderVectors(vector<int>& pasive, vector<int>& orderer)
{
	vector<int> indexes(pasive.size());
	std::iota(indexes.begin(), indexes.end(), 0);
#if _DEBUG
	std::cout << "Ini: Pasive:" << printVector(pasive) << " Ind:" << printVector(indexes) << " Orderer:" << printVector(orderer) << "\n";
#endif
	// Order the indexes vector to map pasive positions.
	std::sort(indexes.begin(), indexes.end(), [pasive](int i1, int i2){return pasive[i1] < pasive[i2]; }); //std::cout << printVector(indexes) << "\n";
	// Order Pasive according to weights in Orderer using indexes mapping
	std::sort(pasive.begin(), pasive.end(), [orderer, indexes](int i1, int i2){return orderer[indexes[i1]] < orderer[indexes[i2]]; }); //std::cout << printVector(pasive) << "\n";
	// simply sort by order Orderer
	std::sort(orderer.begin(), orderer.end());
#if _DEBUG
	std::cout << "End: Pasive:" << printVector(pasive) << " Ind:" << printVector(indexes) << " Orderer:" << printVector(orderer) << "\n";
#endif
}

// Given a pair vector, re-order it in same order as the normal vector indexes
void orderPairVectorWithIds(vector<pair<int, int>>& toOrder, vector<int>& indexes)
{
	vector< pair< int, int > > temp;
	int permSize = indexes.size() - 1;
	for (int j = 0; j <= permSize; j++)
	{
		int job = indexes[j];
		int position = -1;
		for (int k = 0; k <= permSize; k++)
		{
			if (toOrder[k].first == job)
				position = k;
		}

		temp.push_back(pair<int, int>(toOrder[position].first, toOrder[position].second));
	}

	toOrder = temp;
}


// Given a vector with values {1,2,3,2} returns a jagged vector {0}{0,0}{0,0,0}{0,0}
vector< vector < int > > createJaggedVector(vector<int> machinesPerStage)
{
	vector< vector< int > > machineFreeingTimes_S_M;
		//= vector< vector< int > >(machinesPerStage.size());

	for (int i = 0; i < machinesPerStage.size(); i++)
	{
		vector<int> temp = vector<int>(machinesPerStage[i], 0);
		machineFreeingTimes_S_M.push_back(temp);
	}

	return machineFreeingTimes_S_M;
}


/*	Function to compute makespan for a Hybrid Flowshop, right now is just a test that takes 
all data like is a simple FS and considers 2-3 stages with same data (m machines per stage). */
long int PfspInstance::computeHMS(vector<int> &sol)
{
	return PfspInstance::computeHMS(sol, nbJob);
}

/* Compute the Total Completion Time of a given solution */
long int PfspInstance::computeHTCT(vector< int > & sol)
{
	return computeHTCT(sol, nbJob);
}

/* Compute the weighted tardiness of a given solution */
long int PfspInstance::computeHWT(vector< int > & sol)
{
	return computeHWT(sol, nbJob);
}

/* Compute the weighted earliness of a given solution */
long int PfspInstance::computeHWE(vector< int > & sol)
{
	return computeHWE(sol, nbJob);
}

/* Compute the weighted earliness tardines of a given solution */
long int PfspInstance::computeHWET(vector< int > & sol)
{
	return computeHWET(sol, nbJob);
}

/* Compute the weighted earliness tardines with Due Date Windows of a given solution */
long int PfspInstance::computeHWETDDW(vector< int > & sol)
{
	return computeHWETDDW(sol, nbJob);
}

#ifdef HORIZONTAL
/*	Function to compute makespan (PARTIAL SOLUTION) for a Hybrid Flowshop */
long int PfspInstance::computeHMS(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	//// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	//// This JaggedVector will have a ghost stage 0 never used... To keep things like other people.
	//vector< int > machinesLine(nbMac, 0);
	//// Jagged Vector with the finishing times of all machines. 
	//vector< vector< int > > machineFreeingTimes_S_M;
	//machineFreeingTimes_S_M = vector< vector< int > >(stages + 1, machinesLine);
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int PT; // Processing Time
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime; 

	for (j = 1; j <= permSize; ++j)
	{
		jobNumber = sol[j];
		previousTaskFT = 0;

		for (i = 1; i <= stages; i++)
		{
			//// Get First Available Machine and Finishing time.
			//FAM_ID = minValue(machineFreeingTimes_S_M[i]);

			// Get First Available Machine improved to reduce idle times 
			//(If there are some machines that end before job can start, we take the one that gives less idle time)
			FAM_ID = maxLower(machineFreeingTimes_S_M[i], previousTaskFT);

			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID];

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];
			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;
			previousTaskFT = endingTime;
		}
	}

	string x = printJaggedVector(machineFreeingTimes_S_M);// , stages + 1, machinesPerStage);
	// Makespan
	return maxValue(machineFreeingTimes_S_M[stages]);// Makespan
}

/*compute partial weighted tardiness*/
long int PfspInstance::computeHWT(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int PT; // Processing Time
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;

	int wT = 0;
	for (j = 1; j <= permSize; ++j)
	{
		jobNumber = sol[j];
		previousTaskFT = 0;
#ifdef _DEBUG
		string xxx = "J" + std::to_string(j) + "(" + std::to_string(jobNumber) + ") PT:" + printVector(processingTimesMatrix[jobNumber]) + "\n";
		xxx += "S0:" + printJaggedVector(machineFreeingTimes_S_M, stages + 1, machinesPerStage) + "\n";
#endif
		// for all stages but last one.
		for (i = 1; i < stages; i++)
		{
			// Get First Available Machine improved to reduce idle times 
			FAM_ID = maxLower(machineFreeingTimes_S_M[i], previousTaskFT);

			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID];

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];
			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;
			previousTaskFT = endingTime;
		}
		// For last stage (has been acumultated already in for (i++).
		// Get First Available Machine improved to reduce idle times 
		FAM_ID = maxLower(machineFreeingTimes_S_M[i], previousTaskFT);

		previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID];

		// Task can't start until has finished in previous stage and until there is a free machine in actual stage
		startingTime = std::max(previousMachineFT, previousTaskFT);
		endingTime = startingTime + processingTimesMatrix[jobNumber][i];
		// Increase finishing times
		machineFreeingTimes_S_M[i][FAM_ID] = endingTime;
		previousTaskFT = endingTime;
#ifdef _DEBUG
		xxx += "S" + std::to_string(i) + ":\n" + printJaggedVector(machineFreeingTimes_S_M, stages + 1, machinesPerStage) + "\n";
		xxx += "Objective: " + std::to_string(wT) + " + max(" + std::to_string(endingTime) + " - "
			+ std::to_string(dueDates[jobNumber]) + ", 0) * " + std::to_string(priority[jobNumber]);
		printf((xxx + "\n****************************************\n").c_str());
#endif
		wT += (std::max(endingTime - dueDates[jobNumber], 0L) * priority[jobNumber]);
	}
	// Weighted Tardiness
	return wT;
}
#else
// OLD VERSION WITH OPTION TO USE 2 VECTORS (PERMUTATION and FINISHING TIMES)
/*	Function to compute makespan (PARTIAL SOLUTION) for a Hybrid Flowshop 
long int PfspInstance::computeHMS(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

#ifdef TWO_VECTORS
	// With the pair (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));
#else
	// Previous task finishing times. 
	vector< int > PTFT = vector< int >(sol.size(), 0);
	// New deep copy of the permutation, so we can play with it
	vector< int > perm = vector< int >(sol);
#endif

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int PT; // Processing Time
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;

	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		for (j = 1; j <= permSize; ++j)
		{
#ifdef TWO_VECTORS
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time
#else
			jobNumber = perm[j]; // Get Job number
			previousTaskFT = PTFT[j]; // Get previous task finishing time
#endif
			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]);// Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;
#ifdef TWO_VECTORS
			jobAndFT[j].second = endingTime;
#else
			PTFT[j] = endingTime;
#endif
		}
#ifdef TWO_VECTORS
		// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
		std::sort(jobAndFT.begin(), jobAndFT.end(),
			[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
			return left.second < right.second; });
#else
		// Rearrange both vectors (permu and finishing times) according to finishing times values. That way for 
		// next stage the permutation will be different, and will be ordered by realease times  in previous stage.
		orderVectors(perm, PTFT);
#endif		
	}
#ifdef _DEBUG
	string x = printJaggedVector(machineFreeingTimes_S_M);
#endif
	// Makespan
	return maxValue(machineFreeingTimes_S_M[stages]);// Makespan
}*/

inline void sortPairVector(vector< pair< int, int > > &jobAndFT)
{
	std::sort(jobAndFT.begin(), jobAndFT.end(),
		[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
		return left.second < right.second; });
}

inline void sortPairVectorSlackTieBreak(vector< pair< int, int > > &jobAndFT, std::vector< long int > &dueDates)
{

	// Cant make it work. 
	/*
	std::sort(jobAndFT.begin(), jobAndFT.end(),
		[](const std::pair<int, int> &left, const std::pair<int, int> &right, std::vector< long int > &dueDates)
	{
		
		if (left.second == right.second)
			return (dueDates[left.first] - left.second) < (dueDates[right.first] - right.second);
		return left.second < right.second; 
	});*/
}

/*	Function to compute makespan (PARTIAL SOLUTION) for a Hybrid Flowshop */
long int PfspInstance::computeHMS(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	// With the pair (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;

	//string XXX;
	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		for (j = 1; j <= permSize; ++j)
		{
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time

			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]);// Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;

			jobAndFT[j].second = endingTime;
		}
		//XXX += printPairedVector(jobAndFT) + "\n";
		// Sounds crazy but doesnt improve... if (i != stages)
		// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
		std::sort(jobAndFT.begin(), jobAndFT.end(),
			[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
			return left.second < right.second; });
			
	}
#ifdef _DEBUG
	string x = printJaggedVector(machineFreeingTimes_S_M);
#endif
	// Makespan
	return maxValue(machineFreeingTimes_S_M[stages]);// Makespan
}

/*	Function to compute Total Completion Time (PARTIAL SOLUTION) for a Hybrid Flowshop */
long int PfspInstance::computeHTCT(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	// With the pair (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;

	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		for (j = 1; j <= permSize; ++j)
		{
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time

			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]);// Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;

			jobAndFT[j].second = endingTime;
		}
		// Sounds crazy but doesnt improve... if (i != stages)
		// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
		std::sort(jobAndFT.begin(), jobAndFT.end(),
			[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
			return left.second < right.second; });
	}

	int TCT = 0;
	for (std::vector< pair< int, int > >::iterator it = jobAndFT.begin(); it != jobAndFT.end(); ++it)
		TCT += it->second;

#ifdef _DEBUG
	string x = printJaggedVector(machineFreeingTimes_S_M);
#endif
	// Makespan
	return TCT;// Makespan
}

/*compute partial weighted tardiness*/
long int PfspInstance::computeHWT(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	// With the pair (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;
	long int wT = 0;
	bool lastStage = false;
	//string tmp = "";

	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		if (i == stages) lastStage = true;
#ifdef _DEBUG
		string tmp = printPairedVector(jobAndFT) + "\n";
		tmp += "S" + std::to_string(i) + " FT: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		for (j = 1; j <= permSize; ++j)
		{
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time
#ifdef _DEBUG
			tmp += "J" + std::to_string(j) + "(" + std::to_string(jobNumber) + ") PTFT: " + std::to_string(previousTaskFT) + " PT:" + std::to_string(processingTimesMatrix[jobNumber][i]) + "\n";
#endif
			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]);// Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;

			jobAndFT[j].second = endingTime;
			if (lastStage) wT += (std::max(endingTime - dueDates[jobNumber], 0L) * priority[jobNumber]);
#ifdef _DEBUG
			tmp += "S" + std::to_string(i) + " FT: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		}
		if (!lastStage){
			// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
			std::sort(jobAndFT.begin(), jobAndFT.end(),
				[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
				return left.second < right.second; });
		}
	}
	// Weighted Tardiness
	return wT;
}

/*compute partial weighted earliness*/
long int PfspInstance::computeHWE(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	// With the pair (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;
	long int wE = 0;
	bool lastStage = false;
	//string tmp = "";

	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		if (i == stages) lastStage = true;
#ifdef _DEBUG
		string tmp = printPairedVector(jobAndFT) + "\n";
		tmp += "S" + std::to_string(i) + " FT: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		for (j = 1; j <= permSize; ++j)
		{
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time
#ifdef _DEBUG
			tmp += "J" + std::to_string(j) + "(" + std::to_string(jobNumber) + ") PTFT: " + std::to_string(previousTaskFT) + " PT:" + std::to_string(processingTimesMatrix[jobNumber][i]) + "\n";
#endif
			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]);// Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;

			jobAndFT[j].second = endingTime;
			if (lastStage) wE += (std::max(dueDates[jobNumber] - endingTime, 0L) * weightsE[jobNumber]);
#ifdef _DEBUG
			tmp += "S" + std::to_string(i) + " FT: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		}
		if (!lastStage){
			// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
			std::sort(jobAndFT.begin(), jobAndFT.end(),
				[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
				return left.second < right.second; });
		}
	}
	// Weighted Earliness
	return wE;
}

/*compute partial weighted earliness tardiness*/
long int PfspInstance::computeHWET(vector<int> &sol, int size)
{
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	// With the tuple (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;
	long int wET = 0;
	bool lastStage = false;
	//string tmp = "";

	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		if (i == stages) lastStage = true;
#ifdef _DEBUG
		string tmp = printPairedVector(jobAndFT) + "\n";
		tmp += "S" + std::to_string(i) + " FT: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		for (j = 1; j <= permSize; ++j)
		{
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time
#ifdef _DEBUG
			tmp += "J" + std::to_string(j) + "(" + std::to_string(jobNumber) + ") PTFT: " + std::to_string(previousTaskFT) + " PT:" + std::to_string(processingTimesMatrix[jobNumber][i]) + "\n";
#endif
			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]);// Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;

			jobAndFT[j].second = endingTime;
			if (lastStage)
			{
				wET += (std::max(dueDates[jobNumber] - endingTime, 0L) * weightsE[jobNumber]); // Earliness
				wET += (std::max(endingTime - dueDates[jobNumber], 0L) * priority[jobNumber]); // Tardiness
			}
#ifdef _DEBUG
			tmp += "S" + std::to_string(i) + " FT: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		}
		if (!lastStage){
			// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
			std::sort(jobAndFT.begin(), jobAndFT.end(),
				[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
				return left.second < right.second; });
		}
	}
	// Weighted Earliness Tardiness
	return wET;
}

#ifdef IDLE_TWO_STEPS
#else
#endif

long int inline simpleIdleInsertion(vector<vector<pair<int, int>>>& LSCompleteSR, 
	vector< int >& idleTI, vector< int >& machinesPerStage, vector<vector<long int>>& processingTimesMatrix,
	vector<long int>& earlinessDD, vector<long int>& tardinessDD, vector<long int>& weightsE)
{
	int stages = machinesPerStage.size() - 1;
	int mInLastStage = machinesPerStage[stages];

	// Add dummy job behind to make the + inifinite (or similar)
	for (int i = 0; i < mInLastStage; i++)
		LSCompleteSR[i].push_back(pair<int, int>(0, std::numeric_limits<int>::max()));

	int improvements = 0;

	// Foreach machine
	for (int i = 0; i < mInLastStage; i++)
	{
		int dimSize = LSCompleteSR[i].size();
		// Foreach job in this machine, starting from the second to last.
		for (int j = dimSize - 2; j >= 0; j--)
		{
			int k = j + 1; // His follower
			int jobj = LSCompleteSR[i][j].first;
			int jobk = LSCompleteSR[i][k].first;

			// Gap is the distance between finishing one job and starting next one. 
			int startk = (LSCompleteSR[i][k].second - processingTimesMatrix[jobk][stages]) + idleTI[jobk];
			int endj = LSCompleteSR[i][j].second + idleTI[jobj];

			// Distance to tardiness due date, negative means, there is not tardiness yet.
			int tardiness = endj - tardinessDD[jobj];
			int earliness = earlinessDD[jobj] - endj;

			int gap = startk - endj;

			// If there is a gap between jobs, and job j has slack to tardiness due date... we move. 
			if (gap > 0 && tardiness < 0)
			{
				// Instead of earliness I use negative tardiness because Its convinient to move the job as much 
				// as posible in his due date window, to avoid doing it again for the next jobs 
				int idleInsertion = std::min(gap, -1 * tardiness);

				idleTI[jobj] += idleInsertion;

				// The idle time inserted that is minor to earliness, after weigthed, is the reduction of the Objective function.
				improvements += (std::max(std::min(earliness, idleInsertion), 0) * weightsE[jobj]);

			}
		}
	}
	return improvements;
}

/*compute partial weighted earliness tardiness with Due Date Windows*/
long int PfspInstance::computeHWETDDW(vector<int> &sol, int size)
{
	// 613 613
	//sol = { 0, 9, 1, 4, 5, 6, 8, 2, 3, 7, 10 }; // Best for Instance0
	
	//// 906 1009 (906 con hack)
	//sol = { 0, 9, 2, 6, 4, 7, 8, 1, 5, 3, 10 }; // Best for Instance1 
	//vector<int> sol2ndStage = { 0, 9, 2, 4, 6, 7, 1, 8, 3, 10, 5 }; // Best for Instance1
	
	// 1041 1573 (1103 si lo dejamos ir solo: sol : 0, 10, 3, 5, 9, 2, 8, 1, 4, 6, 7 ) (1041 con hack )
	//sol = { 0, 9, 3, 10, 2, 5, 1, 4, 8, 7, 6 }; // Best for Instance2
	//vector<int> sol2ndStage = { 0,  3, 10, 9, 5, 8, 2, 1, 7, 4, 6 }; // Best for Instance2
	
	// 357 357
	//sol = { 0, 4, 9, 7, 1, 6, 2, 10, 5, 3, 8 }; // Best for Instance3
	
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();
	int mInLastStage = machinesPerStage[stages];

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	// With the pair (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));

	// Machine assigned to the job marked by index. 
	vector< int > LSMA(getNbJob()+1);
	// Idle Time Inserted
	vector< int > idleTI(getNbJob() + 1);

	// Last Stage complete representation. A Vector of pairs <job, finishingTime> for each machine. 
	vector< vector < pair< int, int > > > LSCompleteSR(mInLastStage); 
	
	string tmp;
	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;
	long int wE = 0, wT = 0;
	long int wET = 0;
	bool lastStage = false;

	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		if (i == stages) lastStage = true;
		//// To force a change in second stage... Just to compare with Pan
		//if (i == stages) orderPairVectorWithIds(jobAndFT, sol2ndStage);
		
#ifdef _DEBUG2
		tmp = printPairedVector(jobAndFT) + "\n";
		tmp += "S" + std::to_string(i) + " MachinesFTs: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		for (j = 1; j <= permSize; ++j)
		{
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time

			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]); // Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

#ifdef _DEBUG2
			tmp += "J" + std::to_string(j) + "(" + std::to_string(jobNumber) + ") PTFT: " + std::to_string(previousTaskFT) + " FAMFT: " + std::to_string(previousMachineFT) + " PT:" + std::to_string(processingTimesMatrix[jobNumber][i]) + "\n";
#endif
			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;

			jobAndFT[j].second = endingTime;
			if (lastStage)
			{
				wE = (std::max(earlinessDD[jobNumber] - endingTime, 0L) * weightsE[jobNumber]); // Earliness With DDW
				wT = (std::max(endingTime - tardinessDD[jobNumber], 0L) * priority[jobNumber]); // Tardiness With DDW
				wET += wE + wT;

				// Save the id of the machine used to assign
				LSMA[jobNumber] = FAM_ID;

				// SR direct representation for last stage:
				LSCompleteSR[FAM_ID].push_back(pair<int, int>(jobNumber, endingTime));
			}
#ifdef _DEBUG2
			tmp += "S" + std::to_string(i) + " MachinesFTs: " + printVector(machineFreeingTimes_S_M[i]) + " wE: " + std::to_string(wE) + " wT: " + std::to_string(wT) + "\n";
#endif
		}
		if (!lastStage){
			// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
			std::sort(jobAndFT.begin(), jobAndFT.end(),
				[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
				return left.second < right.second; });
		}

	}

	// Add dummy job behind to make the + inifinite (or similar)
	for (int i = 0; i < mInLastStage; i++)
		LSCompleteSR[i].push_back(pair<int, int>(0, std::numeric_limits<int>::max()));

	int improvements = 0;

#ifdef SIMPLE_IDLE_INSERTION312
	// Foreach machine
	for (int i = 0; i < mInLastStage; i++)
	{
		int dimSize = LSCompleteSR[i].size();
		// Foreach job in this machine, starting from the second to last.
		for (int j = dimSize - 2; j >= 0; j--)
		{
			int k = j + 1; // His follower
			int jobj = LSCompleteSR[i][j].first;
			int jobk = LSCompleteSR[i][k].first;

			// Gap is the distance between finishing one job and starting next one. 
			int startk = (LSCompleteSR[i][k].second - processingTimesMatrix[jobk][stages]) + idleTI[jobk];
			int endj = LSCompleteSR[i][j].second + idleTI[jobj];

			// Distance to tardiness due date, negative means, there is not tardiness yet.
			int tardiness = endj - tardinessDD[jobj];
			int earliness = earlinessDD[jobj] - endj;
			
			int gap = startk - endj;

			// If there is a gap between jobs, and job j has slack to tardiness due date... we move. 
			if (gap > 0 && tardiness < 0)
			{
				// Instead of earliness I use negative tardiness because Its convinient to move the job as much 
				// as posible in his due date window, to avoid doing it again for the next jobs 
				int idleInsertion = std::min(gap, -1 * tardiness);

				idleTI[jobj] += idleInsertion;

				// The idle time inserted that is minor to earliness, after weigthed, is the reduction of the Objective function.
				improvements += (std::max(std::min(earliness, idleInsertion), 0) * weightsE[jobj]);

				//break;
			}
		}
	}
	wET -= improvements;
	improvements = 0;
#else

	improvements = simpleIdleInsertion(LSCompleteSR, idleTI, machinesPerStage, processingTimesMatrix, earlinessDD, tardinessDD, weightsE);
		
	wET -= improvements;
	improvements = 0;

#endif
	
	////////////////////////////////////////////////////////////////////////
	//////   Problemo con los indices... Al meter el job final   ///////////
	//////   para tener el mas infinito, ahora tengo un indice   ///////////
	//////   de std::numeric_limits<int>::max(); que se buscara en los jobs, o se metera...  ///////////
	// Intento de solucion, poniendo 0, ya que será el job fantasma de Federico con 0 PT :)
	////////////////////////////////////////////////////////////////////////


	// Try now to move block is needed. There could be the case:
	// 4 Jobs are together. Job 1 and 2 have earliness with acumulated w of 10, job 3 has no tardiness for next 10 positions 
	// with weight 6 and job 4 has tardiness weight 6. Between jobs 4 and 5 there is a gap of 13. 
	// If we move the block (all jobs), we must consider simultaneously all weight and Tardiness/Earliness, also gap with next job.

	// The way of doing it is studying forward job by job, making a block, and cumulating all weights. 
	// If weights show an improvement posible we can move all the block to the right. But we should do it 
	// a maximum number of positions equivalent to the minimum of all the changeing states of that formula
	// (when we fill the gap with next job or one of the jobs changes of state(early, tardy or on time))

	// Foreach machine
	for (int i = 0; i < mInLastStage; i++)
	{
		int dimSize = LSCompleteSR[i].size();

		// Foreach job in this machine, starting from the second to last.
		for (int j = dimSize - 2; j >= 0; j--)
		{
			//int k = j + 1; // His follower
			int jobj = LSCompleteSR[i][j].first;
			
			int gap = 0;
			int minMovementToChange = std::numeric_limits<int>::max();
			int weightsBalance = 0;
			int endOfBlock = 0;
			
			int k = -1;
			// Elaborate the block checking gaps between elements and his following
			for (k = j; k < dimSize - 1; k++)
			{
				gap = 0;
				int jobk = LSCompleteSR[i][k].first;
				int nextJob = LSCompleteSR[i][k + 1].first;

				// Gap is the distance between finishing one job and starting next one. 
				int endk = LSCompleteSR[i][k].second + idleTI[jobk];
				int nextStart = (LSCompleteSR[i][k + 1].second - processingTimesMatrix[nextJob][stages]) + idleTI[nextJob];

				int movementToChange = 0;

				// Distance to tardiness due date, negative means, there is not tardiness yet.
				int tardiness = endk - tardinessDD[jobk];
				int earliness = earlinessDD[jobk] - endk;

				// If early, save distance to end of earliness and Earliness Weight * -1.
				if (earliness > 0)
				{
					movementToChange = earliness;
					weightsBalance += weightsE[jobk] * -1; 
					// Its negative because we are reducing Objective function value when inserting idle times. 
				}
				else
				{
					// Else if tardy. Save distance to tardiness and Tardiness Weight
					if (tardiness >= 0)
					{
						movementToChange = std::numeric_limits<int>::max();
						weightsBalance += priority[jobk];
					}
					// Else, If on time, save distance to tardiness and weight 0.
					else
					{
						movementToChange = tardiness * -1;
						weightsBalance += 0;
					}
				}

				gap = nextStart - endk;

				minMovementToChange = std::min(minMovementToChange, movementToChange);

				if (gap > 0)
				{
					endOfBlock = k;
					minMovementToChange = std::min(minMovementToChange, gap);
					
					// We break for, because we already have a block.
					break;
				}
			}

			// Now we have a block j to k (j =< k)
			// We have the movement until the balance of weights change. 
			// And we have the change of the objective function per unit moved (can be 0 or negative)

			// If movement is not beneficial, continue studing next element. 
			// if (weightsBalance >= 0) continue;
			// If movement  
			if (weightsBalance > 0) continue;
			else
			{
				// Move the elements with idle time insertion
				for (int x = j; x <= k; x++)
				{
					int jobx = LSCompleteSR[i][x].first;
					idleTI[jobx] += minMovementToChange;
				}

				improvements += minMovementToChange * weightsBalance;

				// BECAUSE WE HAD IMPROVED, BUT COULD BE REQUIRED ANOTHER MOVE IN SAME BLOCK TO IMRPOVE MORE
				// I DO THIS DIRTY SHIT UNTILL I SWAP TO DoWhile with conditional incrementation
				j = j++;
			}

		}
	}

	wET += improvements;

	// Weighted Earliness Tardiness
	return wET;
}







/*compute partial weighted earliness tardiness with Due Date Windows*/
/*
long int PfspInstance::computeHWETDDW(vector<int> &sol, int size)
{
	sol = {0, 9, 1, 4, 5, 6, 8, 2, 3, 7, 10};
	int permSize = size;
	// number of stages
	int stages = getNbStages();

	// Stages as a vector with number of machines per stage.
	vector< int > machinesPerStage = getStages();

	// Create a jagged vector [Stages][Machines] to store the finishing times of each machine in each stage. 
	vector< vector< int > > machineFreeingTimes_S_M = createJaggedVector(machinesPerStage);

	// With the pair (jobId, PreviousStageFinishingTimes);
	vector< pair< int, int > > jobAndFT;
	for (int j = 0; j <= permSize; j++)
		jobAndFT.push_back(pair<int, int>(sol[j], 0));

	// Machine assigned to the job marked by index. 
	vector< int > LSMA(size + 1);
	// Idle Time Inserted
	vector< int > idleTI(size + 1);

	int j, i; // indexes job/stage
	int jobNumber; // jobID
	int FAM_ID; // First Available Machine
	int previousTaskFT, previousMachineFT;
	int endingTime, startingTime;
	long int wE = 0, wT = 0;
	long int wET = 0;
	bool lastStage = false;
	//string tmp = "";

	// In first stage all our previous task finishing times are 0. (no previous task)
	// In consecuent stages we will use a pair based vector with (jobId, finishingTime).
	for (i = 1; i <= stages; ++i)
	{
		if (i == stages) lastStage = true;
#ifdef _DEBUG
		string tmp = printPairedVector(jobAndFT) + "\n";
		tmp += "S" + std::to_string(i) + " MachinesFTs: " + printVector(machineFreeingTimes_S_M[i]) + "\n";
#endif
		for (j = 1; j <= permSize; ++j)
		{
			jobNumber = jobAndFT[j].first; // Get Job number
			previousTaskFT = jobAndFT[j].second; // Get previous task finishing time

			// TODO integrate both steps in a single function.
			FAM_ID = minValue(machineFreeingTimes_S_M[i]);// Get First Available Machine 
			previousMachineFT = machineFreeingTimes_S_M[i][FAM_ID]; // Get FAM finishing Time.

#ifdef _DEBUG
			tmp += "J" + std::to_string(j) + "(" + std::to_string(jobNumber) + ") PTFT: " + std::to_string(previousTaskFT) + " FAMFT: " + std::to_string(previousMachineFT) + " PT:" + std::to_string(processingTimesMatrix[jobNumber][i]) + "\n";
#endif
			// Task can't start until has finished in previous stage and until there is a free machine in actual stage
			startingTime = std::max(previousMachineFT, previousTaskFT);
			endingTime = startingTime + processingTimesMatrix[jobNumber][i];

			// Increase finishing times
			machineFreeingTimes_S_M[i][FAM_ID] = endingTime;

			jobAndFT[j].second = endingTime;
			if (lastStage)
			{
				//wET += (std::max(dueDates[jobNumber] - endingTime, 0L) * weightsE[jobNumber]); // Earliness
				//wET += (std::max(endingTime - dueDates[jobNumber], 0L) * priority[jobNumber]); // Tardiness
				wE = (std::max(earlinessDD[jobNumber] - endingTime, 0L) * weightsE[jobNumber]); // Earliness With DDW
				wT = (std::max(endingTime - tardinessDD[jobNumber], 0L) * priority[jobNumber]); // Tardiness With DDW
				wET += wE + wT;

				// Save the id of the machine used to assign
				LSMA[jobNumber] = FAM_ID;
			}
#ifdef _DEBUG
			tmp += "S" + std::to_string(i) + " MachinesFTs: " + printVector(machineFreeingTimes_S_M[i]) + " wE: " + std::to_string(wE) + " wT: " + std::to_string(wT) + "\n";
#endif
		}
		if (!lastStage){
			// Sort pairs Job/FinishingTimes acordint to finishing times, that way we modify permutation for next stage in order of released jobs.
			std::sort(jobAndFT.begin(), jobAndFT.end(),
				[](const std::pair<int, int> &left, const std::pair<int, int> &right) {
				return left.second < right.second; });
		}

	}


	// IDLE time  insertion function  ONLY FOR GAPS...
	int improvements = 0;
	//do
	//{
		improvements = 0;
		// Foreach machine in last stage. 
		for (int m = 0; m < machinesPerStage[stages]; m++)
		{
			for (int j = permSize; j >= 1; j--)
			{
				int gap = 0;
				int jobID = jobAndFT[j].first;

				int machine = LSMA[jobID];
				if (machine != m) continue; // we study one machine each time. Ignore all jobs that are not related to this machine. 

				int endingTime = jobAndFT[j].second + idleTI[jobID];
				int earliness = earlinessDD[jobID] - endingTime;
				int earlinessW = weightsE[jobID];
				int tardiness = endingTime - tardinessDD[jobID];
				int tardinessW = priority[jobID];

				//if (earliness <= 0) continue; // If its not early we dont care about this job anymore. 
				if (tardiness > 0) continue; // If job can be delayd without cost in tardiness... we delay it.

				// Check next job in same machine and move right untill colides with next job or
				// reaches the tardiness due date (even if the job has no tardiness we move it anyway)
				for (int k = j + 1; k <= permSize; k++)
				{
					int kID = jobAndFT[k].first;

					int machinek = LSMA[kID];
					if (machinek != m) continue; // we study one machine each time. Ignore all jobs that are not related to this machine. 

					int endingTimek = jobAndFT[k].second + idleTI[kID];
					int startingTimek = endingTimek - processingTimesMatrix[kID][stages];
					int earlinessk = earlinessDD[kID] - endingTime;
					int earlinessWk = weightsE[kID];
					int tardinessk = endingTime - tardinessDD[kID];
					int tardinessWk = priority[kID];

					if (endingTime < startingTimek)
					{
						gap = startingTimek - endingTime;
					
						//int idleInsertion = std::min(gap, earliness);
						// Instead of earliness I use negative tardiness because Its convinient to move 
						// the job as much as posible in his due date window, to avoid doing it again. 
						int idleInsertion = std::min(gap, -1*tardiness);

						idleTI[jobID] += idleInsertion;

						// The idle time inserted that is minor to earliness, after weigthed, is the reduction of the Objective function.
						improvements += (std::max(std::min(earliness , idleInsertion), 0) * earlinessW);

						break;
					}
				}

			}
		}
	//} while (improvements > 0);

		bool improved;
		//At this point we have the idle times inserted for all gaps... we have to study groups of jobs to see is a idle insertion in x job together improves Earliness more than screws Tardiness.
		do
		{
			improved = false;

			for (int j = 0; j < permSize; j++)
			{

			}


		} while (improved == true);


	// Weighted Earliness Tardiness
	return wET - improvements;
}*/

#endif

/*
long int PfspInstance::computeHWT(vector<int> &sol, int size)
{
	int j;
	long int wt;
	// We need end times on previous machine : 
	vector< long int > previousMachineEndTime(nbJob + 1, 0);

	// Returns the Complition time of all jobs in the vector previousMachineEndTime 
	computePartialMakespans(sol, previousMachineEndTime, processingTimesMatrix, size, nbMac);

	wt = 0;
	for (j = 1; j <= size; ++j)
		wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

	return wt;
}*/

#pragma endregion NEWCODE

long int PfspInstance::computeMS(vector<int> &sol)
{
    int j, m;
    int jobNumber;

    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /* 1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];
    }

    /* others machines : */
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


long int PfspInstance::computeMS(vector<int> &sol,int size)
{
    int j, m;
    int jobNumber;

    // We need end times on previous machine :
    vector< long int > previousMachineEndTime ( nbJob + 1 );
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

/* Compute the weighted tardiness of a given solution starting from a given machine end time table and a starting index */
/**/
long int PfspInstance::computeWT(vector< int > & sol, vector< vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i)
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
    int j,m;

    int jobNumber;
    int end_i = nbJob-1;
   // std::vector< std::vector < int >> head(previousMachineEndTimeMatri);
    long int a_h = 0;
    long int a_t = 0;
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

    k = nbJob-1;



       for ( m = 2; m <= nbMac; ++m )
       {
           int n = nbMac-m+1;
           head[m][1] = head[m-1][1] + processingTimesMatrix[sol[1]][m];
           tail[n][k] = tail[n+1][k] + processingTimesMatrix[sol[k]][n];
       }

      for ( j = 2; j <= end_i; ++j )
        {
            k = nbJob-j;
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



void inline computeTailss(std::vector<int> &sol, int size,std::vector< std::vector< int > > & tail,vector< vector< long> >& processingTimesMatrix,int nbMac)
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

/*
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


long int PfspInstance::computeWT(vector<int> &sol,vector<int>& prevJob,int job,vector<int>& previousMachineEndTime)
{
    int j, m;
    int jobNumber;
    long int wt;

    /* And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /* 1st machine : */
    jobNumber = sol[job];
    previousMachineEndTime[job] = prevJob[1] + processingTimesMatrix[jobNumber][1];
    prevJob[1] = previousMachineEndTime[job];// -> qua iniziare ad aggiornare prevjob[machine 1]
    for ( j = job+1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];
    }


    /* others machines : */
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

//this function sets up the prevJob and previousMachineEndTime vector so that they can be used by the function above
void PfspInstance::computeWTs(vector<int> &sol,vector<int>& prevJob,int job,vector<int>& previousMachineEndTime)
{
    int j, m;
    int jobNumber;

    /* And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /* 1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= job; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1];                           
    }
     prevJob[1]= previousMachineEndTime[job];

    /* others machines : */
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

int PfspInstance::computeIdleTimeCoeff(vector<int>& prevJob, int job)
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
/*Compute flowtime*/
long int PfspInstance::computeFT(vector< int >& sol)
{
    /*TO MODIFY IF WE HAVE TO DEAL WITH RELEASE DATES!!!*/
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += priority[sol[j]];

    return wt;
}

long int PfspInstance::computeFT(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;
    for ( j = 1; j<= size; ++j )
        wt += priority[sol[j]];

    return wt;
}

/* Compute the weighted completion time of a given solution */
long int PfspInstance::computeWCT(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (previousMachineEndTime[j]  * priority[sol[j]]);

    return wt;
}

/*compute partial weighted completion time*/
long int PfspInstance::computeWCT(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (previousMachineEndTime[j]  * priority[sol[j]]);
    }

    return wt;

}

/* total completion time*/
long int PfspInstance::computeTCT(vector< int > &sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;

    for ( j = 1; j<= nbJob; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

long int PfspInstance::computeTCT(vector< int > &sol,int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<=size; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

/* Compute the weighted earliness of a given solution */
long int PfspInstance::computeWE(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);


    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);

    return wt;
}

/*compute partial weighted tardiness*/
long int PfspInstance::computeWE(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);


    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);
    }

    return wt;

}

/* Compute the weighted tardiness of a given solution */
long int PfspInstance::computeT(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L));// * priority[sol[j]]);

    return wt;
}

/*compute partial  tardiness*/
long int PfspInstance::computeT(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) );//* priority[sol[j]]);
    }

    return wt;

}

/* Compute the earliness of a given solution */
long int PfspInstance::computeE(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );

    return wt;
}

/*compute partial earliness */
long int PfspInstance::computeE(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialMakespans(sol, previousMachineEndTime,processingTimesMatrix,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );
    }

    return wt;

}

/*No wait permutation flowshop*/

/* compute partials completition time distances*/

inline void computeNoWaitTimeDistances(vector<int> &sol, int nbMac, int size,vector< vector < long int > >& processingTimesMatrix, vector< int > & completionTimeDistance)
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
        for (int k = 1; k< nbMac ; k++)
        {
            int temp_ctd = processingTimesMatrix[i][k];
            for(int h=k; h < nbMac ; h++ )
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

/* makespan */

long int PfspInstance::computeNWMS(vector<int> &sol)
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

/*partial makespan*/

long int PfspInstance::computeNWMS(vector<int> &sol, int size)
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

/* No wait weighted tardiness*/
long int PfspInstance::computeNWWT(vector< int > &sol)
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

long int PfspInstance::computeNWWT(vector< int > &sol, int size)
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

/*No wait weighted earliness*/
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

/* No wait earliness*/

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

/*No wait tardiness*/

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

/* total completion time*/
long int PfspInstance::computeNWTCT(vector< int > &sol)
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

long int PfspInstance::computeNWTCT(vector< int > &sol,int size)
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

long int PfspInstance::computeNWWCT(vector< int > &sol)
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

long int PfspInstance::computeNWWCT(vector< int > &sol,int size)
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

/*No idle permutation flowshop*/


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

    /*for(int m=1;m<nbMac;++m)
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



long int PfspInstance::computeNIMS(vector<int> &sol, long int nims)
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

long int PfspInstance::computeNIMS(std::vector<int> &sol,int size)
{
   long int nims = processingTimesMatrix[sol[1]][1];;

   std::vector< int > minimumDiff(nbMac,0);
   for(int m=1;m<nbMac;++m)
       minimumDiff[m] = processingTimesMatrix[sol[1]][m+1];


   for(int j=2; j<= size; ++j)
   {

       nims += processingTimesMatrix[sol[j]][1];
       for(int m=1;m<nbMac;++m)
       {
           minimumDiff[m] = std::max(minimumDiff[m]-processingTimesMatrix[sol[j-1]][m],0L) + processingTimesMatrix[sol[j]][m+1];

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

long int PfspInstance::computeSDSTMS(vector<int> &sol)
{
    int j, m;
    int jobNumber;

    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /* 1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1] + setUpTimes[1][sol[j-1]][jobNumber];
    }

    /* others machines : */
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

long int PfspInstance::computeSDSTMS(vector<int> &sol, int size)
{
    int j, m;
    int jobNumber;

    // We need end times on previous machine :
    vector< long int > previousMachineEndTime ( nbJob + 1 );
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

            long stpluspme = setUpTimes[m][sol[j-1]][jobNumber] + previousMachineEndTime[j];
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

inline void computePartialSDSTMakespans( vector< int >& sol, vector< long int >& previousMachineEndTime,vector< vector< long> >& processingTimesMatrix,std::vector< std::vector< std::vector < int > > >& setUpTimes,int nbJob, int nbMac)
{
    int j, m;
    int jobNumber;
   /* And the end time of the previous job, on the same machine : */
    long int previousJobEndTime;

    /* 1st machine : */
    previousMachineEndTime[0] = 0;
    for ( j = 1; j <= nbJob; ++j )
    {
        jobNumber = sol[j];
        previousMachineEndTime[j] = previousMachineEndTime[j-1] + processingTimesMatrix[jobNumber][1] + setUpTimes[1][sol[j-1]][jobNumber];
    }

    /* others machines : */
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

/* Compute the weighted completion time of a given solution */
long int PfspInstance::computeSDSTWCT(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (previousMachineEndTime[j]  * priority[sol[j]]);

    return wt;
}

/*compute partial weighted completion time*/
long int PfspInstance::computeSDSTWCT(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (previousMachineEndTime[j]  * priority[sol[j]]);
    }

    return wt;

}

/* total completion time*/
long int PfspInstance::computeSDSTTCT(vector< int > &sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;

    for ( j = 1; j<= nbJob; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

long int PfspInstance::computeSDSTTCT(vector< int > &sol,int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<=size; ++j ){

        wt += previousMachineEndTime[j];
    }

    return wt;
}

/* Compute the weighted earliness of a given solution */
long int PfspInstance::computeSDSTWE(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);


    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);

    return wt;
}

/*compute partial weighted tardiness*/
long int PfspInstance::computeSDSTWE(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);


    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) * priority[sol[j]]);
    }

    return wt;

}

/* Compute the weighted tardiness of a given solution */
long int PfspInstance::computeSDSTT(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L));// * priority[sol[j]]);

    return wt;
}

/*compute partial  tardiness*/
long int PfspInstance::computeSDSTT(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) );//* priority[sol[j]]);
    }

    return wt;

}

/* Compute the earliness of a given solution */
long int PfspInstance::computeSDSTE(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );

    return wt;
}

/*compute partial earliness */
long int PfspInstance::computeSDSTE(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);

    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(dueDates[sol[j]] - previousMachineEndTime[j] , 0L) );
    }

    return wt;

}

/* Compute the weighted tardiness of a given solution */
long int PfspInstance::computeSDSTWT(vector< int > & sol)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 );
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,nbJob,nbMac);

    wt = 0;
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}

/*compute partial weighted tardiness*/
long int PfspInstance::computeSDSTWT(vector<int> &sol, int size)
{
    int j;
    long int wt;
    /* We need end times on previous machine : */
    vector< long int > previousMachineEndTime ( nbJob + 1 ,0);
    /* And the end time of the previous job, on the same machine : */
    computePartialSDSTMakespans(sol, previousMachineEndTime,processingTimesMatrix,setUpTimes,size,nbMac);


    wt = 0;

    for ( j = 1; j<= size; ++j ){

        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);
    }

    return wt;

}




