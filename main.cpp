#include "pch.h"
#include "vector.h"
#include "matrix.h"

int main(int argc, char** argv) {
    using namespace std;
    
    Matrix<double> m1{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Matrix<double> m2{{10, 0, -8}, {2, -4, 3}, {-5, 7, 6}};
    Vector<double> v1{1, 1};
    Vector<double> v2{2, 2, 2};
    Matrix<double> m3(2, 2, 0);

    m1.shape();
    cout << m1.trace() << endl;
    cout << m1+m2 << endl;

    cout << m3 << endl;

    m3.push_back(v1, 0);
    
    cout << m3 << endl;

    m3.push_back(v2, 1);
    
    cout << m3 << endl;


    m3.push_back(m1, 0);

    cout << m3 << endl;

    m3.pop_back(0);

    cout << m3 << endl;

    m3.pop_back(1);

    cout << m3 << endl;

    m3.clear();

    cout << m3 << endl;
    return 0;
}