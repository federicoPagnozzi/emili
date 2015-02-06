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
#include "pfspinstance.h"


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

int PfspInstance::getNbMac()
{
  return nbMac;
}

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

	strcat(fileNameOK, aux2);
    if(!silence)
    {
	cout << "name : " << fileNameOK << endl;
	cout << "file : " << fileName << endl;
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
            cout << "File " << fileName << " is now open, start to read..." << std::endl;
            cout << "Number of jobs : " << nbJob << std::endl;
            cout << "Number of machines : " << nbMac << std::endl;
            cout << "Memory allowed." << std::endl;
            cout << "Start to read matrix..." << std::endl;
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


/* Compute the weighted tardiness of a given solution */
long int PfspInstance::computeWT(vector< int > & sol)
{
	int j, m;
	int jobNumber;
	long int wt;

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

	wt = 0;
	for ( j = 1; j<= nbJob; ++j )
	    wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

	return wt;
}

/*compute partial weighted tardiness*/
long int PfspInstance::computeWT(vector<int> &sol, int size)
{
    int j, m;
    int jobNumber;
    long int wt;
    // We need end times on previous machine :
    vector< long int > previousMachineEndTime ( size + 1 );
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
/*
long int PfspInstance::computeWT(vector< int > & sol, vector< vector<int > >& previousMachineEndTimeMatrix, int start_i, int end_i)
{
   int m = start_i;
   int j;
   long int wt;
   int jobNumber;
   for(j=start_i;j<end_i;j++)
   {
       jobNumber = sol[j];
       previousMachineEndTimeMatrix[1][j] = previousMachineEndTimeMatrix[1][j-1] + processingTimesMatrix[jobNumber][1];
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
*/

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
