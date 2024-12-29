#include "matrix.hpp"

int main(int argc, char* argv[]) {
    string file_name;
    if (argc > 1){
        file_name = argv[1];
    }else{
        file_name = "input";
    }
    ifstream inputFile(file_name);
    if (!inputFile) {
        cerr << "Input file doesn't exist." << endl;
        return 1;
    }

    string line;
    vector<double> x1, x1e, x2, x2e, x3, x3e, x4, x4e, x5, x5e, y, ye;
    int N_X1, N_X2, N_X3, N_X4, N_X5;

    while (getline(inputFile, line)) {
        if (line[0] == '#') continue; 
        istringstream ss(line);
        if (isdigit(line[0])){
            int N_X1Val, N_X2Val, N_X3Val, N_X4Val, N_X5Val;
            ss >> N_X1Val >> N_X2Val >> N_X3Val >> N_X4Val >> N_X5Val;
            N_X1 = N_X1Val;
            N_X2 = N_X2Val;
            N_X3 = N_X3Val;
            N_X4 = N_X4Val;
            N_X5 = N_X5Val;
            cout << "Order -> ";
            cout << "1st:" << N_X1 << "   " << "2nd:" << N_X2 << "   " << "3rd:" << N_X3 << "   " << "4th:" << N_X4 << "   " << "5th:" << N_X5 << endl;
        }else{
            double x1Val, x2Val, x3Val, x4Val, x5Val, yVal;
            ss >> x1Val >> x2Val >> x3Val >> x4Val >> x5Val >> yVal;

            x1.push_back(x1Val);
            x2.push_back(x2Val);
            x3.push_back(x3Val);
            x4.push_back(x4Val);
            x5.push_back(x5Val);
            y.push_back(yVal);
        }
    }
    inputFile.close();

    if(x1.size() != y.size()){
        cerr << "input file is invalid." << endl;
        return 1;
    }

    int N = (N_X1 + 1) * (N_X2 + 1) * (N_X3 + 1) * (N_X4 + 1) * (N_X5 + 1);
    int M = x1.size();
    const int width = 15;
    cout << "A number of input data -> " << M << endl;
    cout << "Input data list:" << endl;
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;
    cout << setw(width) << "X1" << setw(width) << "X2" << setw(width) << "X3" << setw(width) << "X4" << setw(width) << "X5" << setw(width) << "Y" << endl;
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;
    for(size_t i=0;i<x1.size();i++){
        cout << fixed << setprecision(6) << setw(width) << x1[i] << setw(width) << x2[i] << setw(width) << x3[i] << setw(width) << x4[i] << setw(width) << x5[i] << setw(width) << y[i] << endl;
    }
    cout << "-----------------------------------------------------------------------------------------------------------" << endl;

    Matrix X(M, vector<double>(N, 0.0));
    for(int m=0;m<M;m++){
        int n = 0;
        for(int i=0;i<(N_X1+1);i++){
            for(int j=0;j<(N_X2+1);j++){
                for(int k=0;k<(N_X3+1);k++){
                    for(int l=0;l<(N_X4+1);l++){
                        for(int p=0;p<(N_X5+1);p++){
                            X[m][n] = pow(x1[m],i) * pow(x2[m],j) * pow(x3[m],k) * pow(x4[m],l) * pow(x5[m],p);
                            n++;
                        }
                    }
                }
            }
        }
    }

    Matrix Y(M,{0.0});
    for(size_t i=0;i<M;i++){
        Y[i] = {y[i]};
    }
    // cout << "Matrix calculation has started!" << endl;
    Matrix A = T(X);
    // cout << "Step.1 finished!" << endl;
    // printMatrix(A);
    Matrix B = Product(A,X);
    // cout << "Step.2 finished!" << endl;
    // printMatrix(B);
    Matrix C = Inverse(B);
    // cout << "Step.3 finished!" << endl;
    // printMatrix(C);
    Matrix D = Product(C,A);
    // cout << "Step.4 finished!" << endl;
    // printMatrix(D);
    Matrix Preresult = Product(D,Y);
    // cout << "Step.5 finished!" << endl;
    // printMatrix(Preresult);

    cout << "Fitting is done." << endl;
    cout << "-----------------------------------------------------------------------------------------------------------"<< endl;
    cout << "Fitting result" << endl;
    cout << "-----------------------------------------------------------------------------------------------------------"<< endl;
    int n_ = 0;
    for(int i_=0;i_<(N_X1+1);i_++){
        for(int j_=0;j_<(N_X2+1);j_++){
            for(int k_=0;k_<(N_X3+1);k_++){
                for(int l_=0;l_<(N_X4+1);l_++){
                    for(int p_=0;p_<(N_X5+1);p_++){
                        cout << fixed << setprecision(6) << i_ << "  " << j_ << "  " << k_ << "  " << l_ << "  " << p_ << "  " << "  " << Preresult[n_][0] << endl;
                        n_++;
                    }
                }
            }
        }
    }
    cout << "-----------------------------------------------------------------------------------------------------------"<< endl;
    return 0;
}