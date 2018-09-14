#include <iostream>
#include <GaussJordan.h>
#include <Simplex.h>

using namespace std;

double **ab;

void gaussjordan(int);

void simplex(int);

int main() {
    int n, s;

    cout << "AX=B" << endl;
    begin:
    cout << "Input rank of the matrix A" << endl;
    cin >> n;
    cout << "Input matrix A" << endl;

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

    ch:
    cin.get();
    cout << "1. Gauss-Jordan" << endl;
    cout << "2. Simplex" << endl;
    cin >> s;
    switch (s) {
        case 1:
            gaussjordan(n);
            break;
        case 2:
            simplex(n);
            break;
        default:
            cout << "Choose one of the requested!" << endl;
            cout << "Do you want to change input? y/n (else exit)" << endl;
            char r;
            cin >> r;
            if (r == 'y') {
                delete[] ab;
                ab = nullptr;
                goto begin;
            } else if (r == 'n') goto ch;
    }
    return 0;
}

void gaussjordan(int n) {
    GaussJordan gaussjordan(ab, n);
    cout << gaussjordan.getResult() << endl;
}

void simplex(int n) {
    Simplex simplex(ab, n);
}