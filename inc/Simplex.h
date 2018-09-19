//
// Created by VS0x01 on 12.09.2018.
//

#ifndef LINPROG_SIMPLEX_H
#define LINPROG_SIMPLEX_H

#include <iostream>

class Simplex {
    int n, m, cm; //dimensions, copy m
    double **ab, **cab, **st, **fst; // matrix A and B and copy, simplex matrix, fake simplex matrix
    double *a, *fa; //coefficients and fake coefficients
    bool artbasis = false;
    int *ijk; //j of basis vectors (i, j, k, ...)
    int *ijksti; //order of basis vectors (i, j, k, ...) in st (i)
    int *ijkstj; //order of basis vectors (i, j, k, ...) in st (j)
    int resi, resj, result;  //resolving i, j, results j
    const int values = n; // f valuation

    void simplexTable();

    void fakeSimplexTable(int m);

    bool isBasisVector(int j);

    int findBasisVector(int j);

    double valueCount(int j);

    void setResolving();

    bool checkBasis();

    void fakeBasis();

public:
    Simplex(double **ab, int n, int m);

    double getResult() const;
};


#endif //LINPROG_SIMPLEX_H
