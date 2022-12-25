#include <iostream>
#include <cmath>
#include <cstring>
#include <cblas.h>
#include <chrono>
#include <vector>


using namespace std;


void print_matrix(double *A, int n) {
    cout << "[";
    for(int i = 0 ; i < n; i++) {
        cout << "[";
        for (int j = 0 ; j < n; j++) {
            cout << A[j * n + i] << ", "; 
        }
        cout << "],";
        cout << endl;
    }
    cout << "]";
    cout << endl;
}


void random_init(double *A, int n) {
    for(int i = 0; i < n * n; i++) {
        A[i] = int(rand()) % 100;
    }
}  


double Frob(double *A, double *B, int n, int m) {
    double error = 0;
    for (int i = 0 ; i < n; i++) {
        for (int j = 0 ; j < m; j++) {
            error += (A[j * n + i] - B[j * n + i]) * (A[j * n + i] - B[j * n + i]);     
        }
    }
    return pow(error, 0.5);
}


void matmul_checker(double *A, double *B, double *C, int n) { // C = AB
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A, n, B, n, 0, C, n);
}


double check_unitary(double *A, double *C, int n) {
    memset(C, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
	C[i + i * n] = 1;
    }
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
    		n, n, n, -1.0, A, n, A, n, 1.0, C, n);
    double error = cblas_dnrm2(n * n, C, 1) / pow(n, 0.5);
    return error;
}


void intraortho(double *Q, double *R, int n, int s) {
    double norm_q_1 = cblas_dnrm2(n, Q, 1);
    R[0] = norm_q_1;
    cblas_dscal(n, 1 / norm_q_1, Q, 1);
    for (int i = 1; i < s; i++) {	
	
	cblas_dgemv(CblasColMajor, CblasTrans, n, i, 1.0, Q, n, Q + n * i, 1, 1.0, R + n * i, 1);
	cblas_dgemv(CblasColMajor, CblasNoTrans, n, i, -1.0, Q, n, R + n * i, 1, 1.0, Q + n * i, 1);
        intraortho(Q + n * i, R + n * i + i, n, 1);
    }
}


void QR(double *Q, double *R, int n, int s) { 
    intraortho(Q, R, n, s);
    for (int k = 1; k <= n / s - 1; k++) {
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
                      k * s, s, n, 1.0, Q, n, Q + n * s * k, n, 1.0, R + n * s * k, n);        
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      n, s, k * s, -1.0, Q, n, R + n * s * k, n, 1.0, Q + n * s * k, n);
        intraortho(Q + n * s * k, R + n * s * k + s * k, n, s);
    }
}


int main() {

    vector <int> sizes = {256, 512, 768, 1024, 1280, 1536, 1792, 2048, 2560,
     			  2816, 3072, 3414, 3756, 4096};
    
    for(int i = 0; i < sizes.size(); i++) {
        int n = sizes[i];
        int s = 128;
        cout << n << endl;
        double *A = (double *) malloc(n * n * sizeof(double));
        random_init(A, n);
        double *Q = (double *) malloc(n * n * sizeof(double));
        memcpy(Q, A, n * n * sizeof(double));
        double *R = (double *) malloc(n * n * sizeof(double)); 
    
        auto begin = std::chrono::steady_clock::now();
        QR(Q, R, n, s);
        auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        std::cout << "The time: " << elapsed_ms.count() << " ms\n";
    }
    
    //double *C = (double *) malloc(n * n * sizeof(double));
    //matmul_checker(Q, R, C, n);
    //cout << "||A - QR|| = " << Frob(C, A, n, n) << endl; 
    //cout << "Unitary: " << check_unitary(Q, C,  n) << endl; 
    
    return 0;
}

