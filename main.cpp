#include <iostream>
#include <GaussJordan.h>
#include <Simplex.h>

using namespace std;

double **ab;

void gaussjordan(int);

void simplex(int, int);

int main() {
    int n, m, s;
    cout << "1. Gauss-Jordan" << endl;
    cout << "2. Simplex" << endl;
    cin >> s;
    cout << "AX=B" << endl;
    switch (s) {
        case 1:
            cout << "Input rank of the matrix A" << endl;
            cin >> n;
            cout << "Input matrix A" << endl;
            gaussjordan(n);
            break;
        case 2:
            cout << "Input dimensions (n x m) of the matrix A" << endl;
            cin >> n >> m;
            simplex(n, m);
            break;
        default:
            cout << "Choose one of the requested!" << endl;
    }
    return 0;
}

void gaussjordan(int n) {
    ab = new double *[n]; //rows (i)
    for (int i = 0; i < n; ++i) {  //columns (j)
        ab[i] = new double[n + 1];
        for (int j = 0; j < n; ++j) {
            cin >> ab[i][j];
        }
    }

    cout << "Input matrix B" << endl;
    for (int k = 0; k < n; ++k) {
        cin >> ab[k][n];
    }

    GaussJordan gaussjordan(ab, n);
    cout << gaussjordan.getResult() << endl;
}

void simplex(int n, int m) {
    cout << "Input matrix A" << endl;

    ab = new double *[n]; //rows (i)
    for (int i = 0; i < n; ++i) {  //columns (j)
        ab[i] = new double[m + 1];
        for (int j = 0; j < m; ++j) {
            cin >> ab[i][j];
        }
    }

    cout << "Input matrix B" << endl;
    for (int k = 0; k < n; ++k) {
        cin >> ab[k][m];
    }

    /*std::cout << std::endl;
    for (int i = 0; i < n; ++i) {
        for (int l = 0; l < m; ++l) {
            std::cout << ab[i][l] << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;*/
    Simplex simplex(ab, n, m);
    cout << "f(X) -> max = "<< simplex.getResult() << endl;
}