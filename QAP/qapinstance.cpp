#include "qapinstance.h"


QAPInstance::QAPInstance() {
    n = 0;
}


QAPInstance::QAPInstance(string QAPLibFile) {

    /* temporary attributes */
    int _n = -1;
    vector< vector< matrixEl > > _A;
    vector< vector< matrixEl > > _B;

    int count = 0, i = 0, j = 0;
    int tmp;
    string line;
    ifstream infile(QAPLibFile);

    vector < string > str(1);

    bool is_n_set = false,
         is_a_set = false;


    if (infile) {

        while (getline(infile, line)) {

            if (line.length() == 0) {
                continue;
            } else {
                if (!is_n_set) {
                    // reading n
                    _n = stoi(line);
                    is_n_set = true;
                    _A.resize(_n);
                    _B.resize(_n);
                    str.resize(_n, "0");
                } else if(!is_a_set) {
                    // reading A
                    str[count++] = line;
                    if (count == _n) {
                        for (i = 0 ; i < _n ; i++) {
                            istringstream stream(str[i]);
                            for (j = 0 ; j < _n ; j++) {
                                stream >> tmp;
                                _A[i].push_back(tmp);
                            }
                        }
                        is_a_set = true;
                        count = 0;
                        str.resize(1, "0");
                        str.resize(_n, "0");
                    }
                } else {
                    // reading B
                    str[count++] = line;
                    if (count == _n) {
                        for (i = 0 ; i < _n ; i++) {
                            istringstream stream(str[i]);
                            for (j = 0 ; j < _n ; j++) {
                                stream >> tmp;
                                _B[i].push_back(tmp);
                            }
                        }
                        count = 0;
                    }
                }
            } // end if line.length == 0

        } // end while (getline(...))

    } // end if (infile)

    n = _n;

    //A.resize(n);
    A = _A;
    //B.resize(n);
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
            stra.append(to_string(this->A[i][j]));
            stra.append(" ");
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

    return str;
}

float QAPInstance::computeObjectiveFunction(QAPSolution* solution) {

    int i, j;
    float value = 0.0;
    vector< int > sol = solution->getSolution();

    for (i = 0 ; i < this->n ; i++) {
        for (j = 0 ; j < this->n ; j++) {
            value += this->A[i][j] * this->B[sol[i]][sol[j]];
        }
    }

    return value;

}
