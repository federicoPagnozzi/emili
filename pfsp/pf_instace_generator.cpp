#include <iostream>
#include <random>

using namespace std;
int maino(int argc, char *argv[]) {
	
	if(argc > 3)
	{
		int rndseed = atoi(argv[3]);
		std::mt19937 generator(rndseed);
		std::uniform_int_distribution<int> distr(1,99);	
		int njobs = atoi(argv[1]);
		int nmacs = atoi(argv[2]);
		std::cout << njobs << " " << nmacs << std::endl;
		int max_j = 0;
		int max_i[nmacs];
		
		for(int i=0;i<njobs;i++)
		{
			int temp_j=0;
			for(int j=1;j<=nmacs;j++)
			{				
				int pt = distr(generator);
				std::cout << j << " " << pt << " ";
				max_i[j-1]=i==0?pt:pt+max_i[j-1];
				temp_j+=pt;
			}
			max_j = temp_j>max_j?temp_j:max_j;
			std::cout << std::endl;
		}
		int maxx_i = 0;
		for(int i=0;i<nmacs;i++)
		{
			maxx_i = max_i[i]>maxx_i?max_i[i]:maxx_i;
		}
		int P=max_j>maxx_i?max_j:maxx_i;
		if(argc > 5)
		{
  			std::cout << "Reldue" << std::endl;
			float T = 0.2;
			float R = 0.6;
		
			T = atof(argv[4]);
			R = atof(argv[5]);
				
			int DDlw = (int)((float)P*(1.0f-T-R/2.0f));
			int DDlh = (int)((float)P*(1.0f-T+R/2.0f));
		
			std::uniform_int_distribution<int> dd(DDlw,DDlh);
			std::uniform_int_distribution<int> weight(1,10);
			for(int i=0;i<njobs;i++)
			{
				std::cout << "-1" << " " << dd(generator) << " -1 " << weight(generator) << std::endl;
			}
		}
		std::cout << "LB " << P << std::endl;
		
	}
	else
	 {
		std::cout << "FLOWSHOP INSTANCE GENERATOR" << std::endl;
		std::cout << "Usage:\n\tgenerator #jobs #machines random_seed [ T_factor DD_range ]" <<std::endl;	
	}
    return 0;
}
