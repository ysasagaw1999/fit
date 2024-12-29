#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using Matrix = std::vector<std::vector<double>>;
using namespace std;
void printMatrix(const Matrix& A);
pair<Matrix, Matrix> LUdecomposition(const Matrix& A);
Matrix T(const Matrix& A);
Matrix Cut(const Matrix& A,size_t i,size_t j);
double Det(const Matrix& A);
Matrix Inverse(const Matrix& A);
Matrix Product(const Matrix& A, const Matrix& B);
#endif