#include <iostream>
#include <GaussJordan.h>

using namespace std;

double **ab;

void gaussjordan(int);

int main() {
    int n;

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
    
    gaussjordan(n);
	
	delete[] ab;
    ab = nullptr;
    return 0;
}

void gaussjordan(int n) {
    GaussJordan gaussjordan(ab, n);
    cout << gaussjordan.getResult() << endl;
}
