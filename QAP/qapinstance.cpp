#include "qapinstance.h"

using namespace emili::qap;
QAPInstance::QAPInstance() {
    n = 0;
    bestKnownValue = -1;
}


QAPInstance::QAPInstance(string QAPLibFile) {

    /* temporary attributes */
    int _n = -1;
    float _bkv = -1;
    vector< vector< matrixEl > > _A;
    vector< vector< matrixEl > > _B;

    int count = 0, i = 0, j = 0;
    long tmp;
    string line;
    ifstream infile;
    infile.open(QAPLibFile);

    vector < string > str(1);

    bool is_n_set = false,
         is_a_set = false;


    if (infile) {
        infile >> _n;
        std::cout << _n << std::endl;
        _A.resize(_n, vector< matrixEl >(_n));
        _B.resize(_n, vector< matrixEl >(_n));
        read_line(infile);
        for (i = 0; i < _n; i++) {
            for (j = 0; j < _n; j++) {
                infile >> _A[i][j];
                //_A[i][j] = _A[i][j] * 10;
            }
        }
        for (i = 0; i < _n; i++) {
            for (j = 0; j < _n; j++) {
                infile >> _B[i][j];
                //_B[i][j] = _B[i][j] * 10;
            }
        }
    } else {
        std::cout << QAPLibFile << "file missing or corrupted, cannot open, exiting" << std::endl;
        exit(1);
    } // end if (infile)

    infile.close();

    n = _n;
    bestKnownValue = _bkv;
    A_orig = _A;
    B_orig = _B;

    d_symmetric_flag = check_symmetry(_A);
    f_symmetric_flag = check_symmetry(_B);

    std::cout << d_symmetric_flag << " " << f_symmetric_flag << std::endl;

    null_diagonal_flag = check_null_diagonal(_A);
    /* if one matrix has already null diagonal we need not check the other */
    if (!null_diagonal_flag)
        null_diagonal_flag = check_null_diagonal(_B);

    make_symmetric_flag = XOR(d_symmetric_flag, f_symmetric_flag);

    if ( make_symmetric_flag && null_diagonal_flag ) {
        if ( !d_symmetric_flag )
            _A = make_matrix_symmetric(_A);
        else if (!f_symmetric_flag)
            _B = make_matrix_symmetric(_B);
        else {
            std::cerr << "One matrix should have been symmetric" << std::endl;
            exit(1);
        }
   }

   bool dd_symmetric_flag = check_symmetry(_A);
   bool  ff_symmetric_flag = check_symmetry(_B);

    A = _A;
    B = _B;

} // QAPInstance::QAPInstance(string QAPLibFile)


QAPInstance::~QAPInstance() { }


string QAPInstance::toString(void) {
    string str  = to_string(this->n),
           stra = "",
           strb = "";
    str.append("\n");

    int i, j;

    for (i = 0; i < this->n ; i++) {
        for (j = 0; j < this->n ; j++) {
            //stra.append(to_string(this->A_orig[i][j]));
            stra.append(to_string(this->A[i][j]));
            stra.append(" ");
            //strb.append(to_string(this->B_orig[i][j]));
            strb.append(to_string(this->B[i][j]));
            strb.append(" ");
        }
        stra.append("\n");
        strb.append("\n");
    }

    str.append("\n");
    str.append(stra);
    str.append("\n");
    str.append(strb);
    str.append("\n");

    if (this->bestKnownValue > -1) {
        str.append("Best known value: ");
        str.append(to_string(this->bestKnownValue));
        str.append("\n");
    }

    return str;
}

bool QAPInstance::check_null_diagonal(vector< vector< matrixEl > > matrix) {
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


bool QAPInstance::check_symmetry(vector< vector< matrixEl > > matrix) {
    /* 
      FUNCTION:      check whether the Matrix matrix is symmetric
      INPUT:         pointer to the matrix
      OUTPUT:        TRUE if symmetric, otherwise FALSE
      (SIDE)EFFECTS: none
    */
    long int   i, j;
      
    for ( i = 0 ; i < n - 1 ; i++ ) {
        for ( j = i + 1 ; j < n ; j++ ) {
            if( matrix[i][j] != matrix[j][i] )
                return false;
        }
    }
    return true;
}


vector< vector< matrixEl > > QAPInstance::make_matrix_symmetric(vector< vector< matrixEl > > matrix) {
    long int i, j, help;

    vector< vector< matrixEl > > m2 = matrix;

    for ( i = 0 ; i < n ; i++ ) {
        for ( j = 0 ; j < i ; j++ ) {
            help = matrix[i][j] + matrix[j][i];
            m2[i][j] = help;
            m2[j][i] = help;
        }
    }

    return m2;
}


double QAPInstance::computeObjectiveFunction(QAPSolution* solution) {

    int i, j;
    double value = 0.0;
    vector< int > sol = solution->getSolution();

    for (i = 0 ; i < this->n ; i++) {
        for (j = 0 ; j < this->n ; j++) {
            value += this->A[i][j] * this->B[sol[i]][sol[j]];
        }
    }

    if (make_symmetric_flag)
        return value / 2;

    return value;

}
