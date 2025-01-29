#include "pch.h"
#include "tensor.h"
#include "matrix.h"

int main(int argc, char** argv) {
    using namespace std;
    
    Matrix<double> m1{{1, 3}, {-2, 2}};
    Matrix<double> m2{{0, 2}, {-2, -3}, {1, 1}};
    Matrix<double> m3{{-1, -2, -2, 3}, {2, 1, -1, 4}, {99, 1, 2, 3}};

    Matrix<double> m4{{1, 1}, {-1, -2}};
    Matrix<double> m5{{-1, 3}, {3, 1}};

    Matrix<double> v1{{1, 2, 3}};

    vector<double> v2{1,2,3,4,5,6,7,8,9};

    Matrix<double> m6(3, 3);
    
    cout << m6 << endl;

    for (auto iter=m6.col_begin(0)+1; iter!=m6.col_end(0)-1; ++iter)
        *iter = 1.;

    cout << m6 << endl;

    for (auto iter=m6.row_begin(0)+1; iter!=m6.row_end(0)-1; ++iter)
        *iter = 2.;

    cout << m6 << endl;
    return 0;
}