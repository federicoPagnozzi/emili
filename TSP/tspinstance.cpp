#include "tspinstance.h"

using namespace emili::tsp;
TSPInstance::TSPInstance() {
    n = 0;
    bestKnownValue = -1;
}


TSPInstance::TSPInstance(string TSPLibFile) {

    /* temporary attributes */
    int _n = -1;
    float _bkv = -1;
    vector< vector< matrixEl > > _D;
    float x, y;

    int count = 0, i = 0, j = 0;
    long tmp;
    string line, buf, ewt;
    ifstream infile;
    infile.open(TSPLibFile);

    vector < string > str(1);

    bool is_n_set = false,
         is_a_set = false;


    if (infile) {
        
        infile >> buf;
    while ( buf.compare("NODE_COORD_SECTION") != 0 ) {
	if ( buf.compare("NAME") == 0 ) {
            infile >> buf;
	    infile >> instance_name;
            std::cout << instance_name << std::endl;
	}
	else if ( buf.compare("NAME:") == 0 ) {
	    infile >> instance_name;
	}
	else if ( buf.compare("COMMENT") == 0 ){
            infile >> buf;
            read_line(infile);
	}
        else if ( buf.compare("COMMENT:") == 0 ){
            read_line(infile);
        }
	else if ( buf.compare("TYPE") == 0 ) {
            infile >> buf;
            infile >> instance_type;
	    if( instance_type.compare("TSP") != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	}
        else if ( buf.compare("TYPE:") == 0 ) {
            infile >> instance_type;
            if( instance_type.compare("TSP") != 0 ) {
                fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
                exit(1);
            }
        }
	else if( buf.compare("DIMENSION") == 0 ){
            infile >> buf;
            infile >> n;
	}
	else if ( buf.compare("DIMENSION:") == 0 ) {
            infile >> n;
	}
	else if( buf.compare("DISPLAY_DATA_TYPE") == 0 ){
            infile >> buf;
            read_line(infile);
	}
	else if ( buf.compare("DISPLAY_DATA_TYPE:") == 0 ) {
            read_line(infile);
	}
	else if( buf.compare("EDGE_WEIGHT_TYPE") == 0 ){
            infile >> buf;
	    infile >> ewt;
	    if ( ewt.compare("EUC_2D") == 0 ) {
	    }
	    else if ( ewt.compare("CEIL_2D") == 0 ) {
	    }
	    else if ( ewt.compare("GEO") == 0 ) {
	    }
	    else if ( ewt.compare("ATT") == 0 ) {
	    }
	    else
		std::cerr << "EDGE_WEIGHT_TYPE" <<  ewt << "not implemented" << std::endl;
	}
	else if( buf.compare("EDGE_WEIGHT_TYPE:") == 0 ){
	    /* set pointer to appropriate distance function; has to be one of 
 * 	       EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
           infile >> ewt;
            if ( ewt.compare("EUC_2D") == 0 ) {
            }
            else if ( ewt.compare("CEIL_2D") == 0 ) {
            }
            else if ( ewt.compare("GEO") == 0 ) {
            }
            else if ( ewt.compare("ATT") == 0 ) {
            }
            else
                std::cerr << "EDGE_WEIGHT_TYPE" <<  ewt << "not implemented" << std::endl;
	}
        infile >> buf;
    }


    if( buf.compare("NODE_COORD_SECTION") == 0 ){
	//trace_print("found section contaning the node coordinates\n");
    } else {
	fprintf(stderr,"\n\nSome error ocurred finding start of coordinates from tsp file !!\n");
	exit(1);
    }

    for ( i = 0 ; i < n ; i++ ) {
        infile >> j >> x >> y;
        coords.push_back(tsp_point{x, y});
    }

    } else {
        std::cout << TSPLibFile << "file missing or corrupted, cannot open, exiting" << std::endl;
        exit(1);
    } // end if (infile)

    infile.close();

    _D.resize(n, vector< matrixEl >(n));

    if ( ewt.compare("EUC_2D") == 0 ) {
        for (i = 0 ; i < n ; i++) {
            for (j = 0  ; j < n ; j++) {
                _D[i][j] = round_distance(i,j);
            }
        }
    } else if ( ewt.compare("CEIL_2D") == 0 ) {
        for (i = 0 ; i < n ; i++) {
            for (j = 0  ; j < n ; j++) {
                _D[i][j] = ceil_distance(i,j);
            }
        }
    } else if ( ewt.compare("GEO") == 0 ) {
        for (i = 0 ; i < n ; i++) {
            for (j = 0  ; j < n ; j++) {
                _D[i][j] = geo_distance(i,j);
            }
        }
    } else if ( ewt.compare("ATT") == 0 ) {
        for (i = 0 ; i < n ; i++) {
            for (j = 0  ; j < n ; j++) {
                _D[i][j] = att_distance(i,j);
            }
        }
    }

    D = _D;

    /*for (i = 0 ; i < n ; i++) {
        for (j = 0 ; j < n ; j++) {
           std::cout << D[i][j] << " ";
        }
        std::cout << std::endl;
    }*/

} // TSPInstance::TSPInstance(string TSPLibFile)


TSPInstance::~TSPInstance() { }


string TSPInstance::toString(void) {
    string str  = to_string(this->n),
           stra = "";
    str.append("\n");

    int i, j;

    for (i = 0; i < this->n ; i++) {
        for (j = 0; j < this->n ; j++) {
            //stra.append(to_string(this->A_orig[i][j]));
            stra.append(to_string(this->D[i][j]));
            stra.append(" ");
        }
        stra.append("\n");
    }

    str.append("\n");
    str.append(stra);
    str.append("\n");

    return str;
}

bool TSPInstance::check_null_diagonal(vector< vector< matrixEl > > matrix) {
    /* 
          FUNCTION:      check whether the Matrix matrix has a zero diagonal
          INPUT:         pointer to the matrix
          OUTPUT:        TRUE if null diagonal, otherwise FALSE
          (SIDE)EFFECTS: none
    */
    long int   i;
  
    for ( i = 0 ; i < n ; i++ ) {
        if( matrix[i][i] != 0 ) {
            return false;
       }
    }
    return true;
}


double TSPInstance::computeObjectiveFunction(TSPSolution* solution) {
    int i, j;
    double value = 0.0;
    vector< int > sol = solution->getSolution();

    for (i = 0 ; i < this->n-1 ; i++) {
        value += this->D[sol[i]][sol[i+1]];
    }
    value += this->D[sol[n-1]][sol[0]];

    return value;
}
