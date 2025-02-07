#include "pch.h"
#include "slice.h"
#include "vector.h"
#include "matrix.h"

int main(int argc, char** argv) {
    using namespace std;
    
    Matrix<double> m1{{10, -8}, {2, 3}};
    Matrix<double> m2{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Vector<double> v1{1, 1};
    Vector<double> v2{2, 1, 5, 2, 2};
    Vector<double> v3;
    Matrix<double> m3(2, 2, 0);

    Slice col_slc(0, 2, 1);
    Slice row_slc(0, 2, 1);


    cout << v2 << endl;

    cout << m2 << endl;

    m2(1, col_slc) = v1;

    v3 = m2(1, col_slc);
    cout << m2(row_slc, 1) << endl;

    cout << v2(col_slc) << endl;

    cout << m3 << endl;
    m3.shape();
    m3.push_back(0);

    cout << m3 << endl;
    m3.shape();
    m3.push_back(1);

    cout << m3 << endl;
    m3.shape();
    m3(row_slc, col_slc) = m1;

    cout << m3 << endl;

    return 0;
}