#include "qapneighborhood.h"
using namespace emili::qap;

QAPNeighborhood::~QAPNeighborhood(void) { }
QAPInsertNeighborhood::~QAPInsertNeighborhood(void) { }

/*emili::Solution* QAPNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}*/

emili::Solution* QAPNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}


/***********************************************
 *                                             *
 *                                             *
 *                                             *
 *           QAPInsertNeighborhood             *
 *                                             *
 *                                             *
 *                                             *
 ***********************************************/


emili::Solution* QAPInsertNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}



void QAPInsertNeighborhood::reset(void) {
    start_position = 0;
    end_position = 1;
    sp_iterations = 0;
    ep_iterations = 1;
    first_iter = true;
}


int QAPInsertNeighborhood::size(void) {
    int n = this->getProblemInstance().getInstance()->getn();
    return (n*(n-1)/2);
}


emili::Neighborhood::NeighborhoodIterator QAPInsertNeighborhood::begin(emili::Solution *base) {
    ep_iterations = 1;
    sp_iterations = 0;
    first_iter = true;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}


emili::Solution* QAPInsertNeighborhood::random(emili::Solution *currentSolution) {
    emili::iteration_increment();
    
    std::vector< int > x(((QAPSolution*)currentSolution)->getSolution());

    int a, b;
    
    a = emili::generateRandomNumber() % n;
    do {
        b = emili::generateRandomNumber() % n;
    } while (a == b);

    int sol_i = x[a];
    x.erase(x.begin()+a);
    x.insert(x.begin()+b,sol_i);


    QAPSolution* _cur = new QAPSolution(x);
    double new_value = problem_instance.evaluateSolution(*_cur);
    _cur->setSolutionValue(new_value);
    
    return(_cur);
}


emili::Solution* QAPInsertNeighborhood::computeStep(emili::Solution* value) {
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

    QAPSolution* _value = (QAPSolution *)value;
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
 *          QAPExchangeNeighborhood            *
 *                                             *
 *                                             *
 *                                             *
 ***********************************************/


void QAPExchangeNeighborhood::reset(void) {
    start_position = 0;
    end_position = 1;
    sp_iterations = 0;
    ep_iterations = 1;
    first_iter = true;
}


emili::Neighborhood::NeighborhoodIterator QAPExchangeNeighborhood::begin(emili::Solution *base) {
    ep_iterations = 1;
    sp_iterations = 0;
    //start_position = 0;
    //end_position = 1;
    first_iter = true;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

int QAPExchangeNeighborhood::size(void) {
    return (n*(n-1)/2);
}


void QAPExchangeNeighborhood::reverseLastMove(emili::Solution *step)
{
    std::vector < int >& newsol = ((QAPSolution*)step)->getSolution();
    std::swap(newsol[start_position],newsol[end_position]);
}


emili::Solution* QAPExchangeNeighborhood::random(emili::Solution *currentSolution) {
    emili::iteration_increment();
    
    std::vector< int > x(((QAPSolution*)currentSolution)->getSolution());

    int a, b;
    
    a = emili::generateRandomNumber() % n;
    do {
        b = emili::generateRandomNumber() % n;
    } while (a == b);

    double newcost = currentSolution->getSolutionValue() + computeDelta(a, b, x);
    int c = x[a];
    x[a] = x[b];
    x[b] = c;

    QAPSolution* _cur = new QAPSolution(x);
    // double dddd = problem_instance.evaluateSolution(*_cur);
    _cur->setSolutionValue(newcost);

    return(_cur);
}


double QAPExchangeNeighborhood::computeDelta(int u, int v, vector< int >& x) {
    double delta = 0.0;
    
    if (!symmetric) {
        /**
         * asymmetric case
         */
        for ( int k = 0 ; k < n ; k++ ) {
            if ( (k != u) && (k != v) ) {
                delta += d[k][u] * ( f[x[k]][x[v]] - f[x[k]][x[u]] ) + 
                         d[k][v] * ( f[x[k]][x[u]] - f[x[k]][x[v]] ) + 
                         d[u][k] * ( f[x[v]][x[k]] - f[x[u]][x[k]] ) + 
                         d[v][k] * ( f[x[u]][x[k]] - f[x[v]][x[k]] );
            }    
        }
        delta += d[u][u] * ( f[x[v]][x[v]] - f[x[u]][x[u]] )+
                 d[u][v] * ( f[x[v]][x[u]] - f[x[u]][x[v]] )+
                 d[v][u] * ( f[x[u]][x[v]] - f[x[v]][x[u]] )+ 
                 d[v][v] * ( f[x[u]][x[u]] - f[x[v]][x[v]] );
    } else {
        /**
         * symmetric case
         */
        for ( int k = 0 ; k < n ; k++ ) {
            if ( (k != u) && (k != v) ) {
                delta += ( d[k][u] - d[k][v] ) * ( f[x[k]][x[v]] - f[x[k]][x[u]] );
            }
        }

        if ( !make_symmetric ) {
            delta = delta * 2;
        }
    }

    return delta;
}


emili::Solution* QAPExchangeNeighborhood::computeStep(emili::Solution* value) {
    //printf("before increment\n");
    emili::iteration_increment();

    //printf("%d %d %d %d %d || ", sp_iterations, ep_iterations, start_position, end_position, first_iter);
    //std::cout << start_position << " " << end_position << " " << sp_iterations << " " << ep_iterations << std::endl;

    if(!first_iter && sp_iterations == 0 && ep_iterations == 1)
    {        
        return nullptr;
    }
    if(ep_iterations < n-1){
        ep_iterations++;
    }
    else
    {
        sp_iterations = (sp_iterations + 1) % (n-1);
        ep_iterations = (sp_iterations + 1) % n;
        start_position = (start_position + 1) % (n-1);
        end_position = start_position;
    }
    end_position = (end_position + 1) % n;

    //printf("%d %d swap( %d , %d ) %d\n", sp_iterations, ep_iterations, start_position, end_position, first_iter);

    /** /end_position = (end_position + 1) % n;
    if (end_position == 0) {
        start_position = (start_position + 1) % n;
        end_position = (start_position + 1) % n;
    }

    printf("startposition %d \n", start_position);

    / * if we returned to (0, 0), then we have explored all of the neighborhood * /
    if (!first_iter && start_position == 0 && end_position == 1) {
      //sp_iterations = 0;
      return nullptr;
    }*/

    first_iter = false;

    QAPSolution* _value = (QAPSolution *)value;
    vector< int >& x = _value->getSolution();

    //printf("compute delta\n");
    double newcost = value->getSolutionValue() + computeDelta(start_position, end_position, x);

    //printf("before swap\n");
    int c = x[start_position];
    x[start_position] = x[end_position];
    x[end_position] = c;
    //printf("after swap\n");

    /*_value->setSolutionValue(_value->getSolutionValue() +
                               computeDelta(start_position,
                                            end_position,
                                            x));*/

    _value->setSolutionValue(newcost);

    return _value;
}


/***********************************************
 *                                             *
 *                                             *
 *                                             *
 *          QAPFirst2optNeighborhood           *
 *                                             *
 *                                             *
 *                                             *
 ***********************************************/


QAPFirst2optNeighborhood::~QAPFirst2optNeighborhood(void) { }


emili::Solution* QAPFirst2optNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}


void QAPFirst2optNeighborhood::reset(void) {
    start_position = 0;
    end_position = 0;
    sp_iterations = 1;
    ep_iterations = 1;
}


int QAPFirst2optNeighborhood::size(void) {
    int n = this->getProblemInstance().getInstance()->getn();
    return (n*(n-1)/2);
}


emili::Neighborhood::NeighborhoodIterator QAPFirst2optNeighborhood::begin(emili::Solution *base) {
    ep_iterations = 1;
    sp_iterations = 1;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}


emili::Solution* QAPFirst2optNeighborhood::random(emili::Solution *currentSolution) {
    if (!first_iter)
        return nullptr;

    return computeStep(currentSolution);
}


emili::Solution* QAPFirst2optNeighborhood::first2opt_symmetric(emili::Solution* _value,
                                                               bool make_symmetric_flag) {

    QAPSolution* value = (QAPSolution*)_value;

    bool improvement = true;
    int  i, j, u, v, k;
    double tmp;
    int original_symmetric_factor; /* = 2: original symmetric instance
                                      = 1: original asymmetric instance 
                                    */

    if ( make_symmetric_flag )
        original_symmetric_factor = 1; /* compensation because of not dividing matrix by 2 */
    else
        original_symmetric_factor = 2;

    improvement = true;

    qap::QAP& probinstance = this->getProblemInstance();
    QAPInstance* instance = probinstance.getInstance();
    
    std::vector< int > q = value->getSolution();
    QAPRandomInitialSolution* initsol = new QAPRandomInitialSolution(probinstance);
    std::vector< int > x = initsol->generateRandomPermutation();
    delete initsol;
    double best_found = value->getSolutionValue();
    int n = x.size();

    QAPSolution* newvalue = new QAPSolution(q);
    newvalue->setSolutionValue(best_found);

    vector< vector< matrixEl > > d = instance->getA();
    vector< vector< matrixEl > > f = instance->getB();

    while ( improvement ) {
        improvement = false;
        for ( i = 0 ; i < n ; i++ ) {
            u = x[i];
            for ( j = 0 ; j < n ; j++ ) {
                v = x[j];
                if (u == v)
                    continue;
                tmp = 0;
                for ( k = 0 ; k < n ; k++ ) {
                    if ( (k != u) && (k != v) ) {
                        tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
                    }    
                }

                if (tmp < 0) {
                    improvement = true;
                    best_found += tmp * original_symmetric_factor;
                    /*for (int lll = 0 ; lll < n ; lll++) {
                        std::cout << q[lll] << " ";
                    }
                    std::cout << " " << std::endl;*/
                    int sw = q[u];
                    q[u] = q[v];
                    q[v] = sw;
                    /*for (int lll = 0 ; lll < n ; lll++) {
                        std::cout << q[lll] << " ";
                    }
                    std::cout << " " << std::endl;*/
                    newvalue->setSolution(q);
                    newvalue->setSolutionValue(best_found);
                }
            }
        }
    }

    return newvalue;
}


emili::Solution* QAPFirst2optNeighborhood::first2opt_asymmetric(emili::Solution* _value,
                                                                bool make_symmetric_flag) {

    QAPSolution* value = (QAPSolution*)_value;

    bool improvement = true;
    int  i, j, u, v, k;
    double tmp;

    improvement = true;

    qap::QAP& probinstance = this->getProblemInstance();
    QAPInstance* instance = probinstance.getInstance();
    
    std::vector< int > q = value->getSolution();
    QAPRandomInitialSolution* initsol = new QAPRandomInitialSolution(probinstance);
    std::vector< int > x = initsol->generateRandomPermutation();
    delete initsol;
    double best_found = value->getSolutionValue();
    int n = x.size();

    QAPSolution* newvalue = new QAPSolution(q);
    newvalue->setSolutionValue(best_found);

    vector< vector< matrixEl > > d = instance->getA();
    vector< vector< matrixEl > > f = instance->getB();

    while ( improvement ) {
        improvement = false;
        for ( i = 0 ; i < n ; i++ ) {
            u = x[i];
            for ( j = 0 ; j < n ; j++ ) {
                v = x[j];
                if (u == v)
                    continue;
                tmp = 0;
                for ( k = 0 ; k < n ; k++ ) {
                    if ( (k != u) && (k != v) ) {
                        tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
                               d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
                               d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
                               d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
                    }    
                }
                tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
                       d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
                       d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
                       d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );

                if (tmp < 0) {
                    improvement = true;
                    best_found += tmp;
                    /*for (int lll = 0 ; lll < n ; lll++) {
                        std::cout << q[lll] << " ";
                    }
                    std::cout << " " << std::endl;*/
                    int sw = q[u];
                    q[u] = q[v];
                    q[v] = sw;
                    /*for (int lll = 0 ; lll < n ; lll++) {
                        std::cout << q[lll] << " ";
                    }
                    std::cout << " " << std::endl;*/
                    newvalue->setSolution(q);
                    newvalue->setSolutionValue(best_found);
                }
            }
        }
    }

    return newvalue;
}


emili::Solution* QAPFirst2optNeighborhood::computeStep(emili::Solution* value) {
    QAPInstance* instance = this->getProblemInstance().getInstance();

    if ( instance->is_made_symmetric() ||
        (instance->is_d_symmetric() && instance->is_f_symmetric() && instance->is_null_diagonal())
       ) {
        return first2opt_symmetric(value, instance->is_made_symmetric());
    } else {
        return first2opt_asymmetric(value, instance->is_made_symmetric());
    }
}




/***********************************************
 *                                             *
 *                                             *
 *                                             *
 *           QAPBest2optNeighborhood           *
 *                                             *
 *                                             *
 *                                             *
 ***********************************************/


QAPBest2optNeighborhood::~QAPBest2optNeighborhood(void) { }


emili::Solution* QAPBest2optNeighborhood::step(emili::Solution *currentSolution) {
    return computeStep(currentSolution);
}


void QAPBest2optNeighborhood::reset(void) {
    start_position = 0;
    end_position = 0;
    sp_iterations = 1;
    ep_iterations = 1;
}


int QAPBest2optNeighborhood::size(void) {
    int n = this->getProblemInstance().getInstance()->getn();
    return (n*(n-1)/2);
}


emili::Neighborhood::NeighborhoodIterator QAPBest2optNeighborhood::begin(emili::Solution *base) {
    ep_iterations = 1;
    sp_iterations = 1;
    return emili::Neighborhood::NeighborhoodIterator(this,base);
}

emili::Solution* QAPBest2optNeighborhood::random(emili::Solution *currentSolution) {
    if (!first_iter)
        return nullptr;

    return computeStep(currentSolution);
}


emili::Solution* QAPBest2optNeighborhood::computeStep(emili::Solution* value) {
    QAPInstance* instance = this->getProblemInstance().getInstance();

    if ( instance->is_made_symmetric() ||
        (instance->is_d_symmetric() && instance->is_f_symmetric() && instance->is_null_diagonal())
       ) {
        return best2opt_symmetric(value, instance->is_made_symmetric());
    } else {
        return best2opt_asymmetric(value, instance->is_made_symmetric());
    }
}


emili::Solution* QAPBest2optNeighborhood::best2opt_symmetric(emili::Solution* _value,
                                                              bool make_symmetric_flag) {

    QAPSolution* value = (QAPSolution*)_value;

    bool improvement = true,
         first_it_flag = true; /* first iteration of local search: TRUE */
    int  i, j, u, v, k;
    double tmp,
           max_decrease; /* largest decrease found so far in neighborhood scan */
    int original_symmetric_factor; /* = 2: original symmetric instance
                                      = 1: original asymmetric instance 
                                    */

    if ( make_symmetric_flag )
        original_symmetric_factor = 1; /* compensation because of not dividing matrix by 2 */
    else
        original_symmetric_factor = 2;

    improvement = true;

    qap::QAP& probinstance = this->getProblemInstance();
    QAPInstance* instance = probinstance.getInstance();
    
    std::vector< int > q = value->getSolution();
    QAPRandomInitialSolution* initsol = new QAPRandomInitialSolution(probinstance);
    std::vector< int > x = initsol->generateRandomPermutation();
    delete initsol;
    double best_found = value->getSolutionValue();
    int n = x.size();

    int rchosen = n, schosen = n;  /* memorize which is best move in current iteration */
    int r, s; /* memorize which is best move in previous iteration */

    QAPSolution* newvalue = new QAPSolution(q);
    newvalue->setSolutionValue(best_found);

    vector< vector< matrixEl > > d = instance->getA();
    vector< vector< matrixEl > > f = instance->getB();

    vector< matrixEl > vtemp(n, 0);
    vector< vector< matrixEl > > move_values(n, vtemp);

    while ( improvement ) {
        improvement = false;
        max_decrease = std::numeric_limits<double>::max();

        /* in the first local search iteration the full neighborhood has to be evaluated */

        if (first_it_flag) {
            first_it_flag = false;
            for ( u = 0 ; u < n-1 ; u++ ) {
                for ( v = 0 ; v < n ; v++ ) {
                    tmp = 0;
                    for ( k = 0 ; k < n ; k++ ) {
                        if ( (k != u) && (k != v) ) {
                            tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
                        }    
                    }

                    tmp *= original_symmetric_factor;

                    move_values[u][v] = tmp;
                    if (tmp < max_decrease) {
                        max_decrease = tmp;
                        rchosen = u;
                        schosen = v;
                    }
                }
            }
        } else {
            for ( u = 0 ; u < n-1 ; u++) {
                for ( v = u+1 ; v < n ; v++) {
                    if (u == r || v == s || u == s || v == r) {
                        tmp = 0;
                        for ( k = 0 ; k < n ; k++ ) {
                            if ( (k != u) && (k != v) ) {
                                tmp += ( d[k][u] - d[k][v] ) * ( f[q[k]][q[v]] - f[q[k]][q[u]] );
                            }    
                        }

                        tmp *= original_symmetric_factor;

                        move_values[u][v] = tmp;
                        if (tmp < max_decrease) {
                            max_decrease = tmp;
                            rchosen = u;
                            schosen = v;
                        }
                    } else {
                        /* change derived from prev iteration, u and v differ from rchosen or schosen */
                        tmp = ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
                              ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] );
                        tmp += move_values[u][v];
                        move_values[u][v] = tmp;
                    }
                    if (tmp < max_decrease) {
                        max_decrease = tmp;
                        rchosen = u;
                        schosen = v;
                    }    
                }
            }
        }

        if ( max_decrease < 0 ) {      /* Obj. function value can be improved */
            assert (rchosen < schosen);
            improvement = true;
            best_found += max_decrease;
            /*for (int lll = 0 ; lll < n ; lll++) {
                std::cout << q[lll] << " ";
            }
            std::cout << " " << std::endl;*/
            int sw = q[rchosen];
            q[rchosen] = q[schosen];
            q[schosen] = sw;
            /*for (int lll = 0 ; lll < n ; lll++) {
                std::cout << q[lll] << " ";
            }
            std::cout << " " << std::endl;*/
            newvalue->setSolution(q);
            newvalue->setSolutionValue(best_found);
            r = rchosen; /* memorize previously done move */
            s = schosen; /* memorize previously done move */
        }
    }

    return newvalue;
}


emili::Solution* QAPBest2optNeighborhood::best2opt_asymmetric(emili::Solution* _value,
                                                              bool make_symmetric_flag) {

    QAPSolution* value = (QAPSolution*)_value;

    bool improvement = true,
         first_it_flag = true; /* first iteration of local search: TRUE */
    int  i, j, u, v, k;
    double tmp,
           max_decrease; /* largest decrease found so far in neighborhood scan */

    improvement = true;

    qap::QAP& probinstance = this->getProblemInstance();
    QAPInstance* instance = probinstance.getInstance();
    
    std::vector< int > q = value->getSolution();
    QAPRandomInitialSolution* initsol = new QAPRandomInitialSolution(probinstance);
    std::vector< int > x = initsol->generateRandomPermutation();
    delete initsol;
    double best_found = value->getSolutionValue();
    int n = x.size();

    int rchosen = n, schosen = n;  /* memorize which is best move in current iteration */
    int r = n, s = n; /* memorize which is best move in previous iteration */

    QAPSolution* newvalue = new QAPSolution(q);
    newvalue->setSolutionValue(best_found);

    vector< vector< matrixEl > > d = instance->getA();
    vector< vector< matrixEl > > f = instance->getB();

    vector< matrixEl > vtemp(n, 0);
    vector< vector< matrixEl > > move_values(n, vtemp);

    while ( improvement ) {
        improvement = false;
        max_decrease = std::numeric_limits<double>::max();

        /* in the first local search iteration the full neighborhood has to be evaluated */

        if (first_it_flag) {
            first_it_flag = false;
            for ( u = 0 ; u < n-1 ; u++ ) {
                for ( v = 0 ; v < n ; v++ ) {
                    tmp = 0;
                    for ( k = 0 ; k < n ; k++ ) {
                        if ( (k != u) && (k != v) ) {
                            tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
                                   d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
                                   d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
                                   d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
                        }    
                    }
                    tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
                           d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
                           d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
                           d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
                    move_values[u][v] = tmp;
                    if (tmp < max_decrease) {
                        max_decrease = tmp;
                        rchosen = u;
                        schosen = v;
                    }
                }
            }
        } else {
            for ( u = 0 ; u < n-1 ; u++) {
                for ( v = u+1 ; v < n ; v++) {
                    if (u == r || v == s || u == s || v == r) {
                        tmp = 0;
                        for ( k = 0 ; k < n ; k++ ) {
                            if ( (k != u) && (k != v) ) {
                                tmp += d[k][u] * ( f[q[k]][q[v]] - f[q[k]][q[u]] ) + 
                                       d[k][v] * ( f[q[k]][q[u]] - f[q[k]][q[v]] ) + 
                                       d[u][k] * ( f[q[v]][q[k]] - f[q[u]][q[k]] ) + 
                                       d[v][k] * ( f[q[u]][q[k]] - f[q[v]][q[k]] );
                            }    
                        }
                        tmp += d[u][u] * ( f[q[v]][q[v]] - f[q[u]][q[u]] )+
                               d[u][v] * ( f[q[v]][q[u]] - f[q[u]][q[v]] )+
                               d[v][u] * ( f[q[u]][q[v]] - f[q[v]][q[u]] )+ 
                               d[v][v] * ( f[q[u]][q[u]] - f[q[v]][q[v]] );
                        
                        move_values[u][v] = tmp;
                        if (tmp < max_decrease) {
                            max_decrease = tmp;
                            rchosen = u;
                            schosen = v;
                        }
                    } else { /* change derived from move_values */
                        tmp = ( d[r][u] - d[r][v] + d[s][v] - d[s][u] ) *
                              ( f[q[s]][q[u]] - f[q[s]][q[v]] + f[q[r]][q[v]] - f[q[r]][q[u]] )
                              + ( d[u][r] - d[v][r] + d[v][s] - d[u][s] ) *
                              ( f[q[u]][q[s]] - f[q[v]][q[s]] + f[q[v]][q[r]] - f[q[u]][q[r]] );
                        tmp += move_values[u][v];
                        move_values[u][v] = tmp;
                    }
                    if (tmp < max_decrease) {
                        max_decrease = tmp;
                        rchosen = u;
                        schosen = v;
                    }
                }
            }
        }

        if ( max_decrease < 0 ) {      /* Obj. function value can be improved */
            assert (rchosen < schosen);
            improvement = true;
            best_found += max_decrease;
            /*for (int lll = 0 ; lll < n ; lll++) {
                std::cout << q[lll] << " ";
            }
            std::cout << " " << std::endl;*/
            int sw = q[rchosen];
            q[rchosen] = q[schosen];
            q[schosen] = sw;
            /*for (int lll = 0 ; lll < n ; lll++) {
                std::cout << q[lll] << " ";
            }
            std::cout << " " << std::endl;*/
            newvalue->setSolution(q);
            newvalue->setSolutionValue(best_found);
            r = rchosen; /* memorize previously done move */
            s = schosen; /* memorize previously done move */
        }
    }

    return newvalue;
}
