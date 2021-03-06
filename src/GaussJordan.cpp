//
// Created by VS0x01 on 12.09.2018.
//

#include <GaussJordan.h>
#include <cmath>
#include <sstream>

GaussJordan::GaussJordan(double **ab, int n) : ab(ab), n(n) {
    this->i = 0;
    this->j = 0;

    for (int k = i; k < n; ++k) {
        setResolving();
        for (int l = 0; l < n; ++l) {
            if (ab[l][j] != 0 && l != resolving) {
                for (int m = n; m >= 0; --m) {
                    ab[l][m] -= ab[resolving][m] * ab[l][j];
                }
            }
        }
        i++;
        j++;

        printMatrix();
    }
}

void GaussJordan::setResolving() {
    resolving = i;
    for (int k = i; k < n; ++k) {
        if (fabs(ab[k][j]) < fabs(ab[resolving][j]) && ab[k][j]) resolving = k;
    }

    if (ab[resolving][j] != 1) {
        for (int k = n; k >= j; --k) {
            ab[resolving][k] = ab[resolving][k] ? ab[resolving][k] / ab[resolving][j] : ab[resolving][k];
        }
    }
    if (resolving != i) {
        std::swap(ab[resolving], ab[i]);
        resolving = i;
    }
}

std::string GaussJordan::getResult() {
    std::ostringstream result;
    result << "Result: ";

    for (int k = 0; k < n; ++k) {
        result << ab[k][n] << " ";
    }
    return result.str();
}

void GaussJordan::printMatrix() {
    for (int k = 0; k < n; ++k) {
        for (int l = 0; l < n + 1; ++l) {
            std::cout << ab[k][l] << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
