#include "pch.h"
#include "vector.h"
#include "matrix.h"

int main(int argc, char** argv) {
    using namespace std;
    
    Matrix<double> m1{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Vector<double> v1{2, -4, 1, 3};
    Matrix<double> m2;

    m1.shape();
    cout << m1.trace() << endl;
    return 0;
}