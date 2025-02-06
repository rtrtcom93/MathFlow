#include "pch.h"
#include "slice.h"
#include "vector.h"
#include "matrix.h"

int main(int argc, char** argv) {
    using namespace std;
    
    Matrix<double> m1{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Matrix<double> m2{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Vector<double> v1{1, 1};
    Vector<double> v2{2, 2, 2};
    Matrix<double> m3(2, 2, 0);

    Slice col_slc(1, 3, 1);
    Slice row_slc(1, 3, 1);
    cout << m1 << endl;
    
    for (auto val = m1.row_begin(1)+1; val != m1.row_end(1); ++val) {
        cout << *val;
    }
    cout << endl;

    cout << m1(1, col_slc) << endl;
    cout << m1(row_slc, 1) << endl;
    cout << m1(row_slc, col_slc) << endl;
    return 0;
}