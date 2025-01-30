#include "pch.h"
#include "vector.h"
#include "matrix.h"

int main(int argc, char** argv) {
    using namespace std;
    
    Matrix<double> m1{{1, 3}, {-2, 2}};
    Matrix<double> m2{{0, 2}, {-2, -3}, {1, 1}};
    Matrix<double> m3{{-1, -2, -2}, {2, 1, -1}};

    Matrix<double> m4{{1, 1}, {-1, -2}};
    Matrix<double> m5{{-1, 3}, {3, 1}};

    Matrix<double> v1{{1, 2, 3}};

    // Vector<double> v2{1,2,3,4,5};
    // Vector<double> v3{5,4,3,2,1};

    Matrix<double> m6(3, 3);
    // Matrix<double> m7;

    cout << m2 * m1 << endl;
    cout << m3 * m2 << endl;
    cout << m2 * m3 << endl;
    cout << m4 * m5 << endl;
    cout << m5 * m4 << endl;
    // cin >> m7;
    cout << m3.diag() << endl;
    cout << m3.upp_diag(1) << endl;
    cout << m3.upp_diag(2) << endl;
    cout << m3.low_diag(1) << endl;

    // cout << m7 << endl;

    // cout << v2 << endl;
    // cout << v2.size() << endl;
    return 0;
}