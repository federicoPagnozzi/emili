#include "tspneighborhood.h"
using namespace emili::tsp;

TSPNeighborhood::~TSPNeighborhood(void) { }
TSPInsertNeighborhood::~TSPInsertNeighborhood(void) { }

/*emili::Solution* TSPNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}*/

emili::Solution* TSPNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}


/***********************************************
 *                                             *
 *                                             *
 *                                             *
 *           TSPInsertNeighborhood             *
 *                                             *
 *                                             *
 *                                             *
 ***********************************************/


emili::Solution* TSPInsertNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}



void TSPInsertNeighborhood::reset(void) {
    start_position = 0;
    end_position = 1;
    sp_iterations = 0;
    ep_iterations = 1;
    first_iter = true;
}


int TSPInsertNeighborhood::size(void) {
    int n = this->getProblemInstance().getInstance()->getn();
    return (n*(n-1)/2);
}


emili::Neighborhood::NeighborhoodIterator TSPInsertNeighborhood::begin(emili::Solution *base) {
    ep_iterations = 1;
    sp_iterations = 0;
    first_iter = true;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}


emili::Solution* TSPInsertNeighborhood::random(emili::Solution *currentSolution) {
    emili::iteration_increment();
    
    std::vector< int > x(((TSPSolution*)currentSolution)->getSolution());

    int a, b;
    
    a = emili::generateRandomNumber() % n;
    do {
        b = emili::generateRandomNumber() % n;
    } while (a == b);

    int sol_i = x[a];
    x.erase(x.begin()+a);
    x.insert(x.begin()+b,sol_i);


    TSPSolution* _cur = new TSPSolution(x);
    double new_value = problem_instance.evaluateSolution(*_cur);
    _cur->setSolutionValue(new_value);
    
    return(_cur);
}


emili::Solution* TSPInsertNeighborhood::computeStep(emili::Solution* value) {
    emili::iteration_increment();

    if(!first_iter && sp_iterations == 0 && ep_iterations == 1) {
        return nullptr;
    }
    if(ep_iterations < n-1){
        ep_iterations++;
    } else {
        sp_iterations = (sp_iterations + 1) % (n-1);
        ep_iterations = (sp_iterations + 1) % n;
        start_position = (start_position + 1) % (n-1);
        end_position = start_position;
    }
    end_position = (end_position + 1) % n;

    TSPSolution* _value = (TSPSolution *)value;
    vector< int >& x = _value->getSolution();

    int sol_i = x[start_position];
    x.erase(x.begin()+start_position);
    x.insert(x.begin()+end_position,sol_i);

    double new_value = problem_instance.evaluateSolution(*_value);
    _value->setSolutionValue(new_value);

    //std::cout << _value->getSolutionRepresentation() << std::endl;
    //std::cout << new_value << std::endl;

    return _value;
//    }
}


/***********************************************
 *                                             *
 *                                             *
 *                                             *
 *          TSPExchangeNeighborhood            *
 *                                             *
 *                                             *
 *                                             *
 ***********************************************/


void TSPExchangeNeighborhood::reset(void) {
    start_position = 0;
    end_position = 1;
    sp_iterations = 0;
    ep_iterations = 1;
    first_iter = true;
}


emili::Neighborhood::NeighborhoodIterator TSPExchangeNeighborhood::begin(emili::Solution *base) {
    ep_iterations = 1;
    sp_iterations = 0;
    //start_position = 0;
    //end_position = 1;
    first_iter = true;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

int TSPExchangeNeighborhood::size(void) {
    return (n*(n-1)/2);
}


void TSPExchangeNeighborhood::reverseLastMove(emili::Solution *step)
{
    std::vector < int >& newsol = ((TSPSolution*)step)->getSolution();
    //std::swap(newsol[start_position],newsol[end_position]);
    //double newcost = ((TSPSolution*)step)->getSolutionValue() + computeDelta(start_position, end_position, newsol);
    int stp = min(start_position, end_position);
    int edp = max(start_position, end_position);
    int c;
    for (int i = 0 ; i <= (edp - stp)/2 ; i++) {
       c = newsol[stp+i];
       newsol[stp+i] = newsol[edp - i];
       newsol[edp - i] = c;
    }
    //((TSPSolution*)step)->setSolutionValue(newcost);
}


emili::Solution* TSPExchangeNeighborhood::random(emili::Solution *currentSolution) {
    emili::iteration_increment();
    
    std::vector< int > x(((TSPSolution*)currentSolution)->getSolution());

    int a, b;
    
    a = emili::generateRandomNumber() % n;
    do {
        b = emili::generateRandomNumber() % n;
    } while (a == b);

    double newcost = currentSolution->getSolutionValue() + computeDelta(a, b, x);
    /*int c = x[a];
    x[a] = x[b];
    x[b] = c;*/
    int stp = min(a, b);
    int edp = max(a, b);
    int c;
    for (int i = 0 ; i <= (edp - stp)/2 ; i++) {
       //std::cout << stp+i << " " << edp - i + 1 << endl;
       c = x[stp+i];
       x[stp+i] = x[edp - i];
       x[edp - i] = c;
    }


    TSPSolution* _cur = new TSPSolution(x);
    // double dddd = problem_instance.evaluateSolution(*_cur);
    _cur->setSolutionValue(newcost);

    return(_cur);
}


double TSPExchangeNeighborhood::computeDelta(int u, int v, vector< int >& x) {
    double delta = 0.0;
   
    /*int s1 = (u+1) % n;
    int s2 = (v+1) % n;
   
    std::cout << u << " " << v <<" " <<  d[x[u]][x[v]] << " " << s1 << " " << s2 << " " << d[x[s1]][x[s2]] <<  " " << d[x[u]][x[s1]] << " " << d[x[v]][x[s2]] << std::endl;
    delta = d[x[u]][x[v]] + d[x[s1]][x[s2]] - d[x[u]][x[s1]] - d[x[v]][x[s2]];
    std::cout << "delta " << delta << std::endl;*/

    int a = min(u,v);
    int b = max(u,v);
    int s1, s2;
    if (a == 0 && b == n-1) {
        delta = 0;
    } else {
        s1 = (a-1 +n) % n;
        s2 = (b+1) % n;
        //std::cout << a << " " << b <<" " <<  d[x[s1]][x[b]] << " " << s1 << " " << s2 << " " << d[x[a]][x[s2]] <<  " " << d[x[s1]][x[a]] << " " << d[x[b]][x[s2]] << std::endl;
        delta = d[x[s1]][x[b]] + d[x[a]][x[s2]] - d[x[s1]][x[a]] - d[x[b]][x[s2]];
    }
    //std::cout << "delta " << delta << std::endl;
 
    return delta;
}


emili::Solution* TSPExchangeNeighborhood::computeStep(emili::Solution* value) {
    //printf("before increment\n");
    emili::iteration_increment();

    //printf("%d %d %d %d %d || ", sp_iterations, ep_iterations, start_position, end_position, first_iter);
    //std::cout << start_position << " " << end_position << " " << sp_iterations << " " << ep_iterations << std::endl;

    if(!first_iter && sp_iterations == 0 && ep_iterations == 1)
    {        
        return nullptr;
    }
    if(ep_iterations < n-1){//n-1
        ep_iterations++;
    }
    else
    {
        /**/sp_iterations = (sp_iterations + 1) % (n-1);
        ep_iterations = (sp_iterations + 1) % n;
        start_position = (start_position + 1) % (n-1);
        end_position = start_position;/**/
        /*sp_iterations = (sp_iterations + 1) % (n-2);
        ep_iterations = (sp_iterations + 1) % (n-1);
        start_position = (start_position + 1) % (n-2);
        end_position = start_position;*/
    }
    end_position = (end_position + 1) % n;

    first_iter = false;

    TSPSolution* _value = (TSPSolution *)value;
    vector< int >& x = _value->getSolution();

    //printf("compute delta\n");
    double newcost = value->getSolutionValue() + computeDelta(start_position, end_position, x);

    //printf("before swap\n");
    /*std::cout << "starting ";
    for (int i = 0 ; i < n ; i++) {
      std::cout << x[i] << " ";
    }
    std::cout << "  |  " << value->getSolutionValue() << std::endl;*/
    int stp = min(start_position, end_position);
    int edp = max(start_position, end_position);
    int c;
    //std::cout << stp << " " << edp  << "    " ;
    for (int i = 0 ; i <= (edp - stp)/2 ; i++) {
       //std::cout << stp+i << " " << edp - i + 1 << endl;
       c = x[stp+i];
       x[stp+i] = x[edp - i];
       x[edp - i] = c;
    }
    /*double sv = 0.0;
    for (int i = 0 ; i < n-1 ; i++) {
       sv += d[x[i]][x[i+1]];
       std::cout << x[i] << " "; //(" << d[x[i]][x[i+1]] << ") ";
    }
    sv += d[x[n-1]][x[0]];
    std::cout << x[n-1];// << " (" << d[x[n-1]][x[0]] << ") ";
    / *for (int i = 0 ; i < n ; i++) {
       std::cout << x[i] << " ";
    }*/
    
    //std::cout << "  |  " << sv << " " << newcost << std::endl;/**/
    //printf("after swap\n");

    /*_value->setSolutionValue(_value->getSolutionValue() +
                               computeDelta(start_position,
                                            end_position,
                                            x));*/

    _value->setSolutionValue(newcost);
    //_value->setSolutionValue(sv);

    return _value;
}


