//
// Created by VS0x01 on 12.09.2018.
//

#include <Simplex.h>
#include <cmath>

Simplex::Simplex(double **ab, int n, int m) : ab(ab), n(n), m(m) {
    bool transform;
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

    if (!checkBasis()) {
        fakeBasis();
    }
    simplexTable();
    transform = transformTable();
    if (artbasis && transform) {
        std::swap(this->m, cm);
        std::swap(st, fst);
        std::swap(this->ab, cab);
        result = m - n;
        artbasis = false;

        for (int j = 0, js = 0, *oi = new int; j < cm - n; ++j) {
            if (!isArtBasisVector(j, oi)) {
                for (int i = 0; i < n; ++i) {
                    st[i][js] = fst[i][j];
                }
                js++;
            } else {
                for (int i = 0; i < m; ++i) {
                    if (vstj[i] > *oi) vstj[i]--;
                }
            }
        }
        for (int i = 0; i < n; ++i) {
            st[i][result] = fst[i][cm - n];
        }
        for (int j = 0; j < m - n + 1; ++j) {
            st[values][j] = fixValue(j);
        }
    }
    transformTable();
    printTable();
}

void Simplex::simplexTable() {
    st = new double *[n + 1];
    for (int i = 0; i < n + 1; ++i) {
        st[i] = new double[m + 1 - n];
    }
    vsti = new int[m];
    vstj = new int[m];
    result = m - n;
    for (int j = 0; j < m; ++j) {
        vsti[j] = -1;
        vstj[j] = -1;
    }
    for (int i = 0; i < n; ++i) {
        vsti[ijk[i]] = i;
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
        st[i][result] = ab[i][m];
    }
    for (int j = 0, js = 0; j < m + 1; ++j) {
        if (!isBasisVector(j)) {
            st[values][js] = valueCount(j);
            vstj[j] = js;
            ++js;
        }
    }
}

bool Simplex::transformTable() {
    while (!checkValues()) {
        printTable();
        setResolving();
        changeJordan();
        if (checkUnsolvability()) {
            std::cout << "No optimal solves" << std::endl;
            return false;
        }
    }
    return true;
}

void Simplex::fakeSimplexTable(int m) {
    cab = new double *[n];
    for (int i = 0; i < n; ++i) {
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

    st = new double *[n + 1];
    for (int i = 0; i < n + 1; ++i) {
        st[i] = new double[this->m + 1 - n];
    }
    std::swap(ab, cab);
    std::swap(this->m, cm);
    std::swap(st, fst);
}

bool Simplex::isBasisVector(int j) {
    for (int i = 0; i < n; ++i) {
        if (j == ijk[i]) return true;
    }
    return false;
}

bool Simplex::isArtBasisVector(int j, int *oi) {
    for (int i = m; i < cm; ++i) {
        if (vstj[i] == j) {
            *oi = j;
            return true;
        }
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

double Simplex::fixValue(int j) {
    double value = 0;
    double a = 0;
    double ai[n];
    for (int i = 0; i < m; ++i) {
        if (vstj[i] == j) a = this->a[i];
    }
    for (int si = 0; si < n; ++si) {
        for (int i = 0; i < m; ++i) {
            if (vsti[i] == si) {
                ai[si] = this->a[i];
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        value += st[i][j] * ai[i];
    }
    return value - a;
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
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (vsti[i] == resi && vstj[j] == resj) {
                std::swap(vsti[i], vsti[j]);
                std::swap(vstj[i], vstj[j]);
                return;
            }
        }
    }
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
    for (int j = 0; j < m; ++j) {
        if (a[j] != 0) {
            std::cout << '\t' << j + 1 << ":" << (vsti[j] != (-1) ? st[vsti[j]][result] : 0) << '\t';
        }
    }
    std::cout << ")" << std::endl;
    return st[n][result];
}

void Simplex::printTable() {
    for (int i = 0; i < m; ++i) {
        std::cout << vsti[i] << '\t' << vstj[i] << "\t" << std::endl;
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
