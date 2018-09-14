//
// Created by VS0x01 on 12.09.2018.
//

#ifndef LINPROG_SIMPLEX_H
#define LINPROG_SIMPLEX_H


#include <iostream>

class Simplex {
    double **ab;
    int n;
public:
    Simplex(double **ab, int n);
};


#endif //LINPROG_SIMPLEX_H
