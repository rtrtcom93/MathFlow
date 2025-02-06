#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "pch.h"
#include "matrix.h"
#include "vector.h"


Vector<double> linspace(double x0, double xn, size_t size);
Vector<double> slice(const Vector<double> &vec, int m, int n);

double norm(const Vector<double> &vec, std::string &&type);

Vector<double> TDMA(const Vector<double>& a, const Vector<double>& b, const Vector<double>& c, const Vector<double>& d);

// void replace(vector<double> &vec1, vector<double> &vec2, size_t m);
// void replace(vector<double> &vec1, vector<double> &&vec2, size_t m);

// void gaussJordan(Matrix<double>& matrix);
// void gaussJordanInv(Matrix<double>& matrix);

#endif