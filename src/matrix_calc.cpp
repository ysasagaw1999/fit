#include "matrix.hpp"

void printMatrix(const Matrix& A) {
    for (const auto& row : A) {
        for (const auto& elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }
}

Matrix T(const Matrix& A){
    size_t n = A.size();
    size_t m = A[0].size();
    Matrix B(m, vector<double>(n, 0.0));
    for(size_t i=0;i<n;i++){
        for(size_t j=0;j<m;j++){
            B[j][i] = A[i][j];
        }
    }
    return B;
}

Matrix Cut(const Matrix& A,size_t i,size_t j){
    size_t n = A.size();
    size_t m = A[0].size();
    Matrix Result(n-1,vector<double>(n-1,0.0));
    size_t i__ = 0;
    for(size_t i_=0;i_<(n-1);i_++){
        if(i_ == i) continue;
        size_t j__ = 0;
        for(size_t j_=0;j_<m;j_++){
            if(j_ == j) continue;
            Result[i__][j__] = A[i_][j_];
            j__++;
        }
        i__++;
    }
    return Result;
}

pair<Matrix, Matrix> LUdecomposition(const Matrix& A){
    size_t n = A.size();
    Matrix L(n,vector<double>(n,0.0));
    Matrix U(n,vector<double>(n,0.0));
    for(size_t i=0;i<n;i++){
        for(size_t j=0;j<n;j++){
            U[i][j] = A[i][j];
            for(size_t k=0;k<i;k++){
                U[i][j] -= L[i][k] * U[k][j]; 
            }
        }
        for(size_t j=i;j<n;j++){
            L[j][i] = A[j][i] / U[i][i]; 
            for(size_t k=0;k<i;k++){
                L[j][i] -= L[j][k] * U[k][i] / U[i][i]; 
            }
        }
    }
    for(size_t i=0;i<n;i++) {
        L[i][i] = 1;
    }
    return {L, U};
}

double Det(const Matrix& A){
    size_t n = A.size();
    size_t m = A[0].size();
    if(n!=m){
        cerr << " you cannot calculate Determinant of a " << n << "x" << m << " matrix." << endl;
        return 0.0;
    }
    if(n==2){
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    }else if(n==3){
        return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    }
    auto [L,U] = LUdecomposition(A);
    double result = 1.0;
    for(size_t i=0;i<n;i++){
        result *= U[i][i];
    }
    return result;
} 

Matrix Product(const Matrix& A, const Matrix& B){
    size_t n = A.size();
    size_t l = A[0].size();
    size_t l_ = B.size();
    size_t m = B[0].size();
    if(l != l_){
        cerr << "You cannot multiply " << n << "x" << l << " to " << l_ << "x" << m << endl;
        return {};
    }
    Matrix Result(n,vector<double>(m,0.0));
    for(size_t i=0;i<n;i++){
        for(size_t j=0;j<m;j++){
            for(size_t k=0;k<l;k++){
                Result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return Result;
}

Matrix Inverse(const Matrix& A){ //LU method
    size_t n = A.size();
    size_t m = A[0].size();
    if(n != m){
        cerr << "You cannot calculate the inverse of a non-square matrix." << endl;
        return {};
    }
    double detA = Det(A);
    if(detA == 0){
        cerr << "Matrix is singular and cannot have an inverse." << endl;
        return {};
    }
    auto [L,U] = LUdecomposition(A);
    Matrix InvU(n,vector<double>(n,0.0));
    Matrix InvL(n,vector<double>(n,0.0));
    for(size_t i=0;i<n;i++){
        InvU[i][i] = 1 / U[i][i];
        InvL[i][i] = 1 / L[i][i];
        for(size_t j=0;j<i;j++){
            for(size_t k=i;k<j;k++){
                InvU[i][j] -= U[i][k] * InvU[k][j] / U[i][i];
            }
        }
        for(size_t j=i+1;j<n;j++){
            for(size_t k=i+1;k<j;k++){
                InvL[i][j] -= L[i][k] * InvL[k][j] / L[i][i];
            }
        }
    }
    return Product(InvU,InvL);
}



