//
// Created by VS0x01 on 12.09.2018.
//

#ifndef LINPROG_SIMPLEX_H
#define LINPROG_SIMPLEX_H

#include <iostream>

class Simplex {
    double **ab, **st;
    double *a;
    int n, m;
    int ijk[3] = {-1, -1, -1}; //j of vector
    int resi, resj, values = n, result; //resolving i, j; f valuation; results j

    void simplexTable();

    double valueCount(int j);

    void setResolving();

    bool checkBasis();

    void fakeBasis();

public:
    Simplex(double **ab, int n, int m);

    double getResult() const;
};


#endif //LINPROG_SIMPLEX_H
