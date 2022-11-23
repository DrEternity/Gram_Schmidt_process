#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>
#include <chrono>
#include <cmath>

int n;

using namespace std;

long double scal_pr(vector <long double> &a, vector <long double> &b, int col_a, int col_b) {
    long double sum = 0;
    for (int i = 0 ; i < n; i++) {
        sum += a[col_a * n + i] * b[col_b * n + i];
    }
    return sum;
}


void step_ort(int k, vector <long double> &A,  vector <long double> &Q, vector <long double> &R, vector <long double> &len_col_Q) { // Шаг ортогонализации [1, n] 
    for (int i = 1 ; i <= k - 1; i++) {
        long double tmp = len_col_Q[i - 1];
        if (tmp == 0) {
            continue;
        }
        long double coef = scal_pr(A, Q, k - 1, i - 1) / tmp; 
        R[(k - 1) * n + (i - 1)] = coef;
        for (int j = 0 ; j < n; j++) {
            Q[(k - 1) * n + j] += coef * Q[(i - 1) * n + j];
        }
    }
    R[(k - 1) * n + (k - 1)] = 1;
    for (int j = 0; j < n; j++) {
        Q[(k - 1) * n + j] *= -1;
        Q[(k - 1) * n + j] += A[(k - 1) * n + j];
        len_col_Q[k - 1] += Q[(k - 1) * n + j] * Q[(k - 1) * n + j];
    }
}


int main(int argc, char* argv[]) // n
{
    auto begin = std::chrono::steady_clock::now(); // системное время    
    ifstream fin("test.txt");
    ofstream fout("res_test.txt");
    n = std::stoi(argv[1]);
    cout << "Размер матрицы " << n << endl; 

    vector <long double > A(n * n);
    for (int i = 0 ; i < n; i++) {
        for (int j = 0 ; j < n; j++) {
            fin >> A[j * n + i]; 
        }
    }
    vector <long double > Q(n * n);
    vector <long double > R(n * n);
    vector <long double> len_col_Q(n);

    for (int i = 1 ; i <= n; i++) { 
        step_ort(i, A, Q, R, len_col_Q);
    }
    for (int i = 0; i < len_col_Q.size(); i++) {
        len_col_Q[i] = pow(len_col_Q[i], 0.5);
    }    

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fout << Q[j * n + i] / len_col_Q[j] << " "; 
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fout << R[j * n + i] * len_col_Q[i] << " ";
        }
    }
    

    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    cout << double (elapsed_ms.count()) / 1000 << " "; 
    fout << double (elapsed_ms.count()) / 1000 << " "; 
    
    fout.close();
    fin.close();
    return 0;
}

