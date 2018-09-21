//
// Created by VS0x01 on 12.09.2018.
//

#include <Simplex.h>
#include <cmath>

Simplex::Simplex(double **ab, int n, int m) : ab(ab), n(n), m(m) {
    /*1 -2 1 0 0
    -2 1 0 1 0
    2 1 0 0 1*/

    /*1 -1 1 0 0
    1 1 0 -1 0
    -1 1 0 0 1*/

    /*-1 1 1 0 0
    1 -1 0 1 0
    1 -2 0 0 1*/

    a = new double[m + 1];
    a[m] = 0;
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
    ijkstj = new int[m];
    nbsti = new int[m];
    nbstj = new int[m];
    result = m - n;
    for (int j = 0; j < m; ++j) {
        ijkstj[j] = -1;
        ijksti[j] = -1;
        nbsti[j] = -1;
        nbstj[j] = -1;
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
    for (int j = 0, js = 0; j < m + 1; ++j) {
        if (!isBasisVector(j)) {
            st[values][js] = valueCount(j);
            nbstj[j] = js;
            ++js;
        }
    }
    while (!checkValues()) {
        printTable();
        setResolving();
        if (checkUnsolvability()) {
            std::cout << "No optimal solves" << std::endl;
            return;
        }
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

bool Simplex::checkValues() {
    tmpresult = st[values][result];
    for (int j = 0; j < m - n; ++j) {
        if (st[values][j] < 0) {
            resj = j;
            return false;
        }
    }
    return true;
}

bool Simplex::checkUnsolvability() {
    return st[values][result] < tmpresult;
}

void Simplex::setResolving() {
    for (int j = 0; j < m - n; ++j) {
        if (fabs(st[values][j]) > fabs(st[values][resj]) && st[values][j] < 0) resj = j;
    }
    double tmp = 0;
    for (int i = 0; i < n; ++i) {
        if (st[i][resj] != 0 && st[i][result] / st[i][resj] > 0) {
            tmp = st[i][result] / st[i][resj];
            resi = i;
            break;
        }
    }
    for (int i = resi + 1; i < n; ++i) {
        if (st[i][resj] != 0 && st[i][result] / st[i][resj] > 0 && st[i][result] / st[i][resj] < tmp) {
            tmp = st[i][result] / st[i][resj];
            resi = i;
        }
    }
    ijksti[ijk[resi]] = -1;
    ijkstj[ijk[resi]] = resj;
    nbsti[resj] = resi;
    nbstj[resj] = -1;
    changeJordan();
}

void Simplex::changeJordan() {
    oldst = new double *[n + 1];
    for (int i = 0; i < n + 1; ++i) {
        oldst[i] = new double[m + 1 - n];
        for (int j = 0; j < m + 1 - n; ++j) {
            oldst[i][j] = st[i][j];
        }
    }
    st[resi][resj] = 1 / st[resi][resj];
    for (int i = 0; i < n + 1; ++i) {
        if (i != resi) st[i][resj] *= (-st[resi][resj]);
    }
    for (int j = 0; j < m + 1 - n; ++j) {
        if (j != resj) st[resi][j] *= (st[resi][resj]);
    }
    for (int i = 0; i < n + 1; ++i) {
        if (i == resi) continue;
        for (int j = 0; j < m + 1 - n; ++j) {
            if (j == resj) continue;
            st[i][j] = (oldst[resi][resj] * oldst[i][j] - oldst[resi][j] * oldst[i][resj]) / oldst[resi][resj];
        }
    }
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
    std::cout << "X* (";
    int rj[n]; // result`s j numbers
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            rj[i] = ijksti[j] == i ? j : nbsti[j] == i ? j : -1;
            if(rj[i] != (-1)) break;
        }
    }
    for (int i = 0; i < n; ++i) {
        if(a[rj[i]] != 0) std::cout << " " << rj[i] + 1 << ":" << st[i][result] << " ";

    }
    std::cout << ")" << std::endl;
    return st[n][result];
}

void Simplex::printTable() {
    for (int i = 0; i < m; ++i) {
        std::cout << ijksti[i] << '\t' << ijkstj[i] << '\t' << nbsti[i] << '\t' << nbstj[i] << std::endl;
    }
    std::cout << std::endl;

    for (int i = 0; i < n + 1; ++i) {
        for (int j = 0; j < m + 1 - n; ++j) {
            std::cout << st[i][j] << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
