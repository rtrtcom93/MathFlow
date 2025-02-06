#include "pch.h"
#include "slice.h"
#include "vector.h"
#include "matrix.h"

int main(int argc, char** argv) {
    using namespace std;
    
    Matrix<double> m1{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Matrix<double> m2{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Vector<double> v1{1, 1};
    Vector<double> v2{2, 2, 2, 2, 2};
    Matrix<double> m3(2, 2, 0);

    Slice col_slc(1, 3, 1);
    Slice row_slc(1, 3, 1);


    cout << v2 << endl;

    cout << m2 << endl;

    m2(1, col_slc) = v1;

    cout << m2 << endl;

    return 0;
}