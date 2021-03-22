#ifndef TSPINSTANCE_H
#define TSPINSTANCE_H


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "tspsolution.h"

namespace emili {
namespace tsp {

using namespace std;

#define XOR(x,y) ((x && !y) || (!x && y))

//#define M_PI 3.14159265358979323846264

typedef long matrixEl;
typedef struct _tsp_point {
    float x;
    float y;
} tsp_point;

/**
 * TSPInstance
 */
class TSPInstance {

private:
    int n;

    string instance_name, instance_type;
    /* distance matrix */
    vector< vector< matrixEl > > D;

    vector< tsp_point > coords;

    vector< int > optPerm;

    TSPSolution *solution;

    float optValue;

    /* bestKnownValue is the value read in Thomas' QAP instances */
    float bestKnownValue;

    /* symmetry checking and handling */
    bool check_null_diagonal(vector< vector< matrixEl > > matrix);

    static double dtrunc (double x) {
        int k;

        k = (int) x;
        x = (double) k;
        return x;
    }

long int  (*distance)(long int, long int);  /* function pointer */

/*    
 *          FUNCTION: the following four functions implement different ways of 
 *                    computing distances for TSPLIB instances
 *          INPUT:    two node indices
 *          OUTPUT:   distance between the two nodes
 * */

long int round_distance (long int i, long int j) 
/*    
 *          FUNCTION: compute Euclidean distances between two nodes rounded to next 
 *                    integer for TSPLIB instances
 *          INPUT:    two node indices
 *          OUTPUT:   distance between the two nodes
 *          COMMENTS: for the definition of how to compute this distance see TSPLIB
 * */
{
    double xd = coords[i].x - coords[j].x;
    double yd = coords[i].y - coords[j].y;
    double r  = sqrt(xd*xd + yd*yd) + 0.5;

    return (long int) r;
}

long int ceil_distance (long int i, long int j) 
/*    
 *          FUNCTION: compute ceiling distance between two nodes rounded to next 
 *                    integer for TSPLIB instances
 *          INPUT:    two node indices
 *          OUTPUT:   distance between the two nodes
 *          COMMENTS: for the definition of how to compute this distance see TSPLIB
 * */
{
    double xd = coords[i].x - coords[j].x;
    double yd = coords[i].y - coords[j].y;
    double r  = sqrt(xd*xd + yd*yd);

    return (long int)(ceil (r));
}

long int geo_distance (int i, int j) 
/*    
 *          FUNCTION: compute geometric distance between two nodes rounded to next 
 *                    integer for TSPLIB instances
 *          INPUT:    two node indices
 *          OUTPUT:   distance between the two nodes
 *          COMMENTS: adapted from concorde code
 *                    for the definition of how to compute this distance see TSPLIB
 * */
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    long int dd;
    float x1 = coords[i].x, x2 = coords[j].x, 
	  y1 = coords[i].y, y2 = coords[j].y;

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (y1);
    min = y1 - deg;
    longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (y2);
    min = y2 - deg;
    longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * std::acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;

}

long int att_distance (long int i, long int j) 
/*    
 *          FUNCTION: compute ATT distance between two nodes rounded to next 
 *                    integer for TSPLIB instances
 *          INPUT:    two node indices
 *          OUTPUT:   distance between the two nodes
 *          COMMENTS: for the definition of how to compute this distance see TSPLIB
 * */
{
    double xd = coords[i].x - coords[j].x;
    double yd = coords[i].y - coords[j].y;
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    long int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}

public:
    TSPInstance(void);
    /**
     * initialize an instance from a TSPLib instance file
     *
     * Please provide an existing && valid file
     */
    TSPInstance(string TSPLibFile);
    ~TSPInstance(void);

    /* Read/Write private attributes */
    int getn(void) {
        return n;
    }

    void setn(int _n) {
        n = _n;
    }

    vector< vector< matrixEl > > getDistanceMatrix(void) {
        return D;
    }

    void setDistanceMatrix(vector< vector< matrixEl > > _D) {
        D = _D;
    }

    TSPSolution* getSolution(void) {
        return solution;
    }

    void setSolution(TSPSolution* _solution) {
        solution = _solution;
    }

    string toString(void);

    double computeObjectiveFunction(TSPSolution* _solution);

    void read_line(ifstream& s) {
        char tmp[10000];
        s.getline(tmp, sizeof(tmp));
    }

}; // TSPInstance

}
}
#endif
