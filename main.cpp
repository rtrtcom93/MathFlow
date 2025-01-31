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

    Matrix<double> v1{1, 2, 3};

    Vector<double> v2{1,2,3,4,5};
    Vector<double> v3{5,4,3,2,1};

    Matrix<double> m6(3, 3);
    // Matrix<double> m7;

    cout<< v2 << endl;

    return 0;
}