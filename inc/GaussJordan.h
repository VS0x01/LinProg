//
// Created by VS0x01 on 12.09.2018.
//

#ifndef LINPROG_GAUSSJORDAN_H
#define LINPROG_GAUSSJORDAN_H

#include <iostream>

class GaussJordan {
    double **ab;
    int n, i, j, resolving;

    void setResolving();

    void printMatrix(); //for faster debugging without gdb

public:
    GaussJordan(double **ab, int n);

    std::string getResult();
};


#endif //LINPROG_GAUSSJORDAN_H
