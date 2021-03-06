//
// Created by VS0x01 on 12.09.2018.
//

#ifndef LINPROG_SIMPLEX_H
#define LINPROG_SIMPLEX_H

#include <iostream>

class Simplex {
    int n, m, cm; //dimensions, copy m
    double **ab, **cab, **st, **fst, **oldst; // matrix A and B and copy, simplex matrix, fake and old simplex matrix
    double *a, *fa; //coefficients and fake coefficients
    double fm = 0; //free member
    bool artbasis = false;
    int *ijk; //j of basis vectors (i, j, k, ...)
    int *vsti; //order of vectors in st (i)
    int *vstj; // -||- (j)
    int resi, resj, result;  //resolving i, j, results j
    double tmpresult = 0;
    const int values = n; // f valuation

    void simplexTable();

    bool transformTable();

    void fakeSimplexTable(int m);

    bool isBasisVector(int j);

    bool isArtBasisVector(int j, int *oi);

    int findBasisVector(int j);

    double valueCount(int j);

    double fixValue(int j);

    bool checkValues();

    bool checkUnsolvability();

    void setResolving();

    void changeJordan();

    bool checkBasis();

    void fakeBasis();

    void printTable();

public:
    Simplex(double **ab, int n, int m);

    double getResult() const;
};


#endif //LINPROG_SIMPLEX_H
