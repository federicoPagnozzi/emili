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
    for ( j = 1; j<= nbJob; ++j )
        wt += (std::max(previousMachineEndTime[j] - dueDates[sol[j]], 0L) * priority[sol[j]]);

    return wt;
}

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




