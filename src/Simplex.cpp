//
// Created by VS0x01 on 12.09.2018.
//

#include <Simplex.h>

Simplex::Simplex(double **ab, int n, int m) : ab(ab), n(n), m(m) {
    /*1 -2 1 0 0
    -2 1 0 1 0
    2 1 0 0 1*/

    /*1 -1 1 0 0
    1 1 0 -1 0
    -1 1 0 0 1*/
    a = new double[m];
    std::cout << "Input vars` coefficients of the objective function f(X) -> max" << std::endl;
    for (int i = 0; i < m; ++i) {
        std::cin >> a[i];
    }
    ijk = new int[n];
    for (int i = 0; i < n; ++i) {
        ijk[i] = -1;
    }

    if (checkBasis()) {
        simplexTable();
    } else {
        fakeBasis();
    }
}

void Simplex::simplexTable() {
    st = new double *[n + 1];
    for (int i = 0; i < n + 1; ++i) {
        st[i] = new double[m + 1 - n];
    }
    ijksti = new int[m];
    ijkstj = new int[m - n];
    for (int i = 0; i < m + 1 - n; ++i) {
        ijkstj[i] = -1;
    }
    for (int i = 0; i < n; ++i) {
        ijksti[ijk[i]] = i;
    }

    for (int j = 0, js = 0; j < m; ++j) {
        if (!isBasisVector(j)) {
            for (int i = 0; i < n; ++i) {
                st[i][js] = ab[i][j];
            }
            js++;
        }
    }
    for (int i = 0; i < n; ++i) {
        st[i][m - n] = ab[i][m];
    }
    for (int j = 0, js = 0; js < m + 1 - n; ++j, ++js) {
        if (!isBasisVector(j)) st[values][js] = valueCount(j);
        else st[values][js] = valueCount(++j);
    }
}

void Simplex::fakeSimplexTable(int m) {
    cab = new double *[n];
    for (int i = 0; i < n; ++i) {  //columns (j)
        cab[i] = new double[m + 1];
        for (int j = 0; j < this->m; ++j) {
            cab[i][j] = ab[i][j];
        }
    }
    for (int j = this->m; j < m; ++j) {
        int ib = findBasisVector(j);
        for (int i = 0; i < n; ++i) {
            cab[i][j] = i == ib ? 1 : 0;
        }
    }
    for (int k = 0; k < n; ++k) {
        cab[k][m] = ab[k][this->m];
    }

    fst = new double *[n];
    for (int i = 0; i < n; ++i) {
        fst[i] = new double[m];
    }
    std::swap(ab, cab);
    std::swap(this->m, cm);
    std::swap(st, fst);
    simplexTable();
}

bool Simplex::isBasisVector(int j) {
    for (int i = 0; i < n; ++i) {
        if (j == ijk[i]) return true;
    }
    return false;
}

int Simplex::findBasisVector(int j) {
    for (int i = 0; i < n; ++i) {
        if (j == ijk[i]) return i;
    }
}

double Simplex::valueCount(int j) {
    double value = 0;
    if (!artbasis) {
        for (int i = 0; i < n; ++i) {
            value += ab[i][j] * a[ijk[i]];
        }
        return value - a[j];
    } else {
        for (int i = 0; i < n; ++i) {
            value += ab[i][j] * fa[ijk[i]];
        }
        return value;
    }
}

void Simplex::setResolving() {

}

bool Simplex::checkBasis() {
    int c = 0;
    int ijktmp = 0;
    double sum;
    for (int j = 0; j < m; ++j) {
        sum = 0;
        for (int i = 0; i < n; ++i) {
            if (ab[i][j] != 0 && ab[i][j] != 1) {
                break;
            } else {
                if (ab[i][j] == 1) ijktmp = i;
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
    artbasis = true;
    std::cout << "F(X, w) -> max" << std::endl;
    int c = 0;
    for (int i = 0; i < n; ++i) {
        if (ijk[i] == (-1)) {
            ijk[i] = m + c;
            c++;
        }
    }
    fa = new double[m + c];
    for (int i = 0; i < m; ++i) {
        fa[i] = 0;
    }
    for (int i = m; i < m + c; ++i) {
        fa[i] = -1;
    }
    cm = m + c;
    fakeSimplexTable(cm);
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
            std::cout << st[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    return 0.1;
}
