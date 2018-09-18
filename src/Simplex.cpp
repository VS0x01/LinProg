//
// Created by VS0x01 on 12.09.2018.
//

#include <Simplex.h>

Simplex::Simplex(double **ab, int n, int m) : ab(ab), n(n), m(m) {
    st = new double *[n + 1];
    for (int i = 0; i < n + 1; ++i) {
        st[i] = new double[m + 1 - n];
    }

    /*1 -2 1 0 0
    -2 1 0 1 0
    2 1 0 0 1*/
    a = new double[m];
    std::cout << "Input vars` coefficients of the objective function f(X) -> max" << std::endl;
    for (int i = 0; i < m; ++i) {
        std::cin >> a[i];
    }
    if (checkBasis()) {
        simplexTable();
    } else {
        fakeBasis();
    }
}

void Simplex::simplexTable() {
    for (int j = 0, js = 0; j < m; ++j) {
        if (j != ijk[0] && j != ijk[1] && j != ijk[2]) {
            for (int i = 0; i < n; ++i) {
                st[i][js] = ab[i][j];
            }
            js++;
        }
    }
    for (int i = 0; i < n; ++i) {
        st[i][m - n] = ab[i][m];
    }
    for (int j = 0; j < m + 1 - n; ++j) {
        st[values][j] = valueCount(j);
    }
}

double Simplex::valueCount(int j) {
    double value = 0;
    for (int i = 0; i < n; ++i) {
        value += ab[i][j]*a[ijk[i]];
    }
    return value-a[j];
}

void Simplex::setResolving() {

}

bool Simplex::checkBasis() {
    short c = 0;
    short ijktmp;
    double sum;
    for (int j = 0; j < m; ++j) {
        sum = 0;
        for (int i = 0; i < n; ++i) {
            if (ab[i][j] != 0 && ab[i][j] != 1) {
                break;
            } else {
                if(ab[i][j] == 1) ijktmp = i;
                sum += ab[i][j];
                if (sum > 1) {
                    break;
                }
            }
            if (i == n - 1) {
                ijk[ijktmp] = j;
                c++;
            }
        }
        if (c == n) return true;
    }
    return false;
}

void Simplex::fakeBasis() {

}

double Simplex::getResult() const {
    /*std::cout << "X*(";
    for (int i = 0; i < n; ++i) {
        std::cout << st[n][result] << ",";

    }
    std::cout << ")" << std::endl;
    return st[n][result];*/
    for (int i = 0; i < n + 1; ++i) {
        for (int j = 0; j < m + 1 - n; ++j) {
            std::cout << st[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return 0.1;
}
